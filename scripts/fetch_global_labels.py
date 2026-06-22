#!/usr/bin/env python3
"""Construit les labels globaux (per-protéine) pour les 1390 protéines ATLAS.

Sortie : Datasets/ATLAS/global_labels/atlas_global_labels.tsv
Colonnes : pdb_chain, uniprot_id, fold_label, tm_label, localization_class,
           disorder_global, acc_mean

Sources :
  - fold_label        : API ATLAS metadata → SCOP_ID[0][0] (classe a/b/c/d)
  - tm_label          : UniProt REST → présence de features "Transmembrane" (binaire)
  - localization_class: UniProt REST → SUBCELLULAR LOCATION, regroupé en 6 classes
  - disorder_global   : mean(RMSF) depuis labels_md/{id}_prod_R1_merged.tsv
  - acc_mean          : mean(acc) depuis data/{id}_dssp.tsv

Note SCOP : le champ JSON `SCOP_class` est la chaîne descriptive
("Alpha and beta proteins (a/b)"). La lettre de classe canonique a/b/c/d est
le 1er token de `SCOP_ID` (ex "c.86.1.1" -> "c"). On utilise donc SCOP_ID.

aggregation_score est volontairement absent : pas de source sans Aggrescan3D
local (mentionné comme limite dans le rapport).

Le script est ré-exécutable : les réponses ATLAS et UniProt sont mises en cache
sur disque pour éviter de re-télécharger à chaque run.
"""
import json
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ID_TRAIN = ROOT / "deciphering" / "id_train.txt"
ID_TEST = ROOT / "deciphering" / "id_test.txt"
LABELS_MD = ROOT / "Datasets" / "ATLAS" / "labels_md"
DSSP_DIR = ROOT / "Datasets" / "ATLAS" / "data"
OUT_DIR = ROOT / "Datasets" / "ATLAS" / "global_labels"
CACHE_DIR = OUT_DIR / "cache"
ATLAS_CACHE = CACHE_DIR / "atlas_metadata.json"
UNIPROT_CACHE = CACHE_DIR / "uniprot.json"
OUT_TSV = OUT_DIR / "atlas_global_labels.tsv"

ATLAS_URL = "https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/metadata/{}"
UNIPROT_BATCH_URL = (
    "https://rest.uniprot.org/uniprotkb/accessions"
    "?accessions={}&fields=accession,ft_transmem,cc_subcellular_location"
)
UNIPROT_CHUNK = 100          # accessions par requête batch
SLEEP_UNIPROT = 1.0          # 1 req/sec max (consigne)
SLEEP_ATLAS = 0.25           # politesse envers l'API ATLAS

# Espèce : on ne garde que les organismes ATLAS >= MIN_SPECIES protéines
# (les autres -> "NA", droppés à l'entraînement, comme fold e/f/g).
# Seuil retenu : 30 prot -> 4 classes (H. sapiens, E. coli, S. cerevisiae,
# T. thermophilus). Choix documenté comme limite : ATLAS est centré sur
# quelques organismes modèles.
MIN_SPECIES = 30

# Regroupement localisation UniProt -> 6 classes (cf. rapport Sofia p.9-10)
# L'ordre de test compte : on prend la 1re localisation qui matche.
LOC_RULES = [
    ("membrane", ["membrane"]),
    ("nucleus", ["nucleus", "nuclear"]),
    ("mitochondrion", ["mitochond"]),
    ("extracellular", ["secreted", "extracellular"]),
    ("cytoplasm", ["cytoplasm", "cytosol"]),
]


def load_ids():
    ids = []
    for f in (ID_TRAIN, ID_TEST):
        for line in f.read_text().splitlines():
            line = line.strip()
            if line:
                ids.append(line)
    # dédoublonnage en gardant l'ordre
    return list(dict.fromkeys(ids))


def load_cache(path):
    if path.exists():
        return json.loads(path.read_text())
    return {}


def save_cache(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data))


def http_get_json(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": "atlas-global-labels/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


# ---------------------------------------------------------------------------
# 1. ATLAS metadata -> SCOP class + UniProt id
# ---------------------------------------------------------------------------
def fetch_atlas(ids, cache):
    missing = [i for i in ids if i not in cache]
    print(f"[ATLAS] {len(missing)} à récupérer ({len(ids) - len(missing)} en cache)")
    for n, pid in enumerate(missing, 1):
        try:
            data = http_get_json(ATLAS_URL.format(pid))
            cache[pid] = data.get(pid, {})
        except urllib.error.HTTPError as e:
            print(f"  [WARN] ATLAS {pid}: HTTP {e.code}")
            cache[pid] = {}
        except Exception as e:
            print(f"  [WARN] ATLAS {pid}: {e}")
            cache[pid] = {}
        if n % 50 == 0:
            print(f"  ATLAS {n}/{len(missing)}")
            save_cache(ATLAS_CACHE, cache)
        time.sleep(SLEEP_ATLAS)
    save_cache(ATLAS_CACHE, cache)
    return cache


def species_keep_set(atlas, ids):
    """Organismes ATLAS avec >= MIN_SPECIES protéines (sur l'ensemble des ids)."""
    from collections import Counter
    c = Counter()
    for pid in ids:
        o = atlas.get(pid, {}).get("organism")
        if o and o != "-":
            c[o] += 1
    keep = {o for o, n in c.items() if n >= MIN_SPECIES}
    print(f"[species] {len(keep)} classes >= {MIN_SPECIES} prot : "
          f"{sorted(keep)}")
    return keep


def species_label(meta, keep):
    o = meta.get("organism")
    if o and o != "-" and o in keep:
        return o
    return None


def scop_class(meta):
    """1er caractère du 1er SCOP_ID -> a/b/c/d. NaN sinon."""
    sid = meta.get("SCOP_ID")
    if isinstance(sid, list) and sid:
        first = str(sid[0]).strip()
        if first and first[0] in "abcdefg":
            return first[0]
    return None


# ---------------------------------------------------------------------------
# 2. UniProt -> tm_label + localization_class
# ---------------------------------------------------------------------------
def fetch_uniprot(accessions, cache):
    todo = [a for a in accessions if a and a not in cache]
    todo = list(dict.fromkeys(todo))
    print(f"[UniProt] {len(todo)} accessions à récupérer ({len(cache)} en cache)")
    for i in range(0, len(todo), UNIPROT_CHUNK):
        chunk = todo[i:i + UNIPROT_CHUNK]
        url = UNIPROT_BATCH_URL.format(",".join(chunk))
        try:
            data = http_get_json(url)
            results = data.get("results", [])
            for rec in results:
                acc = rec.get("primaryAccession")
                if acc:
                    cache[acc] = parse_uniprot(rec)
            # accessions sans résultat (obsolètes) -> marquer vide
            returned = {r.get("primaryAccession") for r in results}
            for a in chunk:
                if a not in returned and a not in cache:
                    cache[a] = {"tm": 0, "loc": None}
        except Exception as e:
            print(f"  [WARN] UniProt chunk {i}: {e}")
            for a in chunk:
                cache.setdefault(a, {"tm": 0, "loc": None})
        print(f"  UniProt {min(i + UNIPROT_CHUNK, len(todo))}/{len(todo)}")
        save_cache(UNIPROT_CACHE, cache)
        time.sleep(SLEEP_UNIPROT)
    return cache


def parse_uniprot(rec):
    # tm_label : au moins une feature Transmembrane
    tm = 0
    for f in rec.get("features", []):
        if f.get("type") == "Transmembrane":
            tm = 1
            break
    # localization : 1re SUBCELLULAR LOCATION
    loc = None
    for c in rec.get("comments", []):
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for sl in c.get("subcellularLocations", []):
                val = sl.get("location", {}).get("value")
                if val:
                    loc = val
                    break
        if loc:
            break
    return {"tm": tm, "loc": loc}


def map_localization(loc):
    if not loc:
        return None
    low = loc.lower()
    for cls, keys in LOC_RULES:
        if any(k in low for k in keys):
            return cls
    return "other"


# ---------------------------------------------------------------------------
# 3. disorder_global (mean RMSF) + acc_mean (mean acc DSSP)
# ---------------------------------------------------------------------------
def mean_rmsf(pid):
    f = LABELS_MD / f"{pid}_prod_R1_merged.tsv"
    if not f.exists():
        return None
    vals = []
    with f.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            idx = header.index("RMSF")
        except ValueError:
            return None
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > idx:
                try:
                    vals.append(float(cols[idx]))
                except ValueError:
                    pass
    return sum(vals) / len(vals) if vals else None


def mean_acc(pid):
    f = DSSP_DIR / f"{pid}_dssp.tsv"
    if not f.exists():
        return None
    vals = []
    with f.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            idx = header.index("acc")
        except ValueError:
            return None
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) > idx:
                try:
                    vals.append(float(cols[idx]))
                except ValueError:
                    pass
    return sum(vals) / len(vals) if vals else None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    ids = load_ids()
    print(f"{len(ids)} protéines à traiter")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    atlas = fetch_atlas(ids, load_cache(ATLAS_CACHE))

    # accessions UniProt depuis metadata ATLAS
    uni_ids = {}
    for pid in ids:
        up = atlas.get(pid, {}).get("UniProt")
        if isinstance(up, str) and up not in ("-", ""):
            uni_ids[pid] = up
    uniprot = fetch_uniprot(list(uni_ids.values()), load_cache(UNIPROT_CACHE))

    keep_species = species_keep_set(atlas, ids)

    rows = []
    stats = {"fold": 0, "tm": 0, "loc": 0, "disorder": 0, "acc": 0, "species": 0}
    for pid in ids:
        meta = atlas.get(pid, {})
        up = uni_ids.get(pid)
        urec = uniprot.get(up, {}) if up else {}

        fold = scop_class(meta)
        tm = urec.get("tm") if up else None
        loc = map_localization(urec.get("loc")) if up else None
        dis = mean_rmsf(pid)
        acc = mean_acc(pid)
        spc = species_label(meta, keep_species)

        if fold is not None: stats["fold"] += 1
        if tm is not None: stats["tm"] += 1
        if loc is not None: stats["loc"] += 1
        if dis is not None: stats["disorder"] += 1
        if acc is not None: stats["acc"] += 1
        if spc is not None: stats["species"] += 1

        rows.append({
            "pdb_chain": pid,
            "uniprot_id": up if up else "NA",
            "fold_label": fold if fold is not None else "NA",
            "tm_label": tm if tm is not None else "NA",
            "localization_class": loc if loc is not None else "NA",
            "disorder_global": f"{dis:.4f}" if dis is not None else "NA",
            "acc_mean": f"{acc:.4f}" if acc is not None else "NA",
            "species_label": spc if spc is not None else "NA",
        })

    cols = ["pdb_chain", "uniprot_id", "fold_label", "tm_label",
            "localization_class", "disorder_global", "acc_mean", "species_label"]
    with OUT_TSV.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"\n[OK] {len(rows)} lignes -> {OUT_TSV}")
    print("Couverture (non-NA) :")
    for k, v in stats.items():
        print(f"  {k:10s}: {v}/{len(ids)} ({100*v/len(ids):.1f}%)")


if __name__ == "__main__":
    main()
