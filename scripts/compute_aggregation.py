#!/usr/bin/env python3
"""Score d'agrégation par protéine — proxy structural type A3D / TANGO (Python 3).

Aggrescan3D (Python 2.7) n'est pas installable sur Apple Silicon. On calcule à la
place un proxy reproductible reposant sur le même principe structural qu'A3D :
l'agrégation est portée par les **patchs hydrophobes / à forte propension β
EXPOSÉS au solvant** (les résidus enfouis ne contribuent pas en l'état replié).

Score par résidu :
    a3d_proxy_i = ( z_KD[aa] + z_Pβ[aa] ) * RSA_i
  - z_KD  : hydrophobicité Kyte-Doolittle, standardisée sur les 20 aa
  - z_Pβ  : propension feuillet β (Chou-Fasman), standardisée
  - RSA_i : accessibilité relative au solvant = acc_DSSP / ASA_max[aa]  (clippée [0,1])
Convention : score élevé = plus propice à l'agrégation.

Score global = moyenne des scores par résidu.

Entrée : Datasets/ATLAS/data/{pdb_chain}_dssp.tsv  (colonnes seq, acc, DSSP_8, DSSP_3)
Cache   : Datasets/ATLAS/aggregation/{pdb_chain}_a3d.csv  (réexécutable)
Sortie  : Datasets/ATLAS/aggregation/aggregation_scores.tsv  (pdb_chain, aggregation_score)

Usage : python compute_aggregation.py [N]      # N = limite (test sur sous-ensemble)
"""
import sys
import csv
from pathlib import Path
from statistics import mean, pstdev

ROOT = Path(__file__).resolve().parent.parent
DSSP_DIR = ROOT / "Datasets" / "ATLAS" / "data"
ID_TRAIN = ROOT / "deciphering" / "id_train.txt"
ID_TEST = ROOT / "deciphering" / "id_test.txt"
OUT_DIR = ROOT / "Datasets" / "ATLAS" / "aggregation"
OUT_TSV = OUT_DIR / "aggregation_scores.tsv"

# Kyte-Doolittle hydropathy (élevé = hydrophobe)
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
      "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
      "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
      "Y": -1.3, "V": 4.2}
# Chou-Fasman β-sheet propensity
PB = {"A": 0.83, "R": 0.93, "N": 0.89, "D": 0.54, "C": 1.19, "Q": 1.10,
      "E": 0.37, "G": 0.75, "H": 0.87, "I": 1.60, "L": 1.30, "K": 0.74,
      "M": 1.05, "F": 1.38, "P": 0.55, "S": 0.75, "T": 1.19, "W": 1.37,
      "Y": 1.47, "V": 1.70}
# ASA maximale théorique (Tien et al. 2013, Å²)
ASA_MAX = {"A": 129, "R": 274, "N": 195, "D": 193, "C": 167, "Q": 225,
           "E": 223, "G": 104, "H": 224, "I": 197, "L": 201, "K": 236,
           "M": 224, "F": 240, "P": 159, "S": 155, "T": 172, "W": 285,
           "Y": 263, "V": 174}


def zscores(scale):
    vals = list(scale.values())
    m, s = mean(vals), pstdev(vals)
    return {k: (v - m) / s for k, v in scale.items()}


Z_KD = zscores(KD)
Z_PB = zscores(PB)


def load_ids():
    ids = []
    for f in (ID_TRAIN, ID_TEST):
        ids += [l.strip() for l in f.read_text().splitlines() if l.strip()]
    return list(dict.fromkeys(ids))


def residue_scores(dssp_path):
    """Retourne la liste des scores proxy par résidu, ou None si illisible."""
    scores = []
    with dssp_path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_seq, i_acc = header.index("seq"), header.index("acc")
        except ValueError:
            return None
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) <= max(i_seq, i_acc):
                continue
            aa = c[i_seq].strip().upper()
            if aa not in KD:
                continue
            try:
                acc = float(c[i_acc])
            except ValueError:
                continue
            rsa = min(max(acc / ASA_MAX[aa], 0.0), 1.0)
            scores.append((Z_KD[aa] + Z_PB[aa]) * rsa)
    return scores


def main():
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else None
    ids = load_ids()
    if limit:
        ids = ids[:limit]
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    rows, ok, fail = [], 0, 0
    fails = []
    for pid in ids:
        dssp = DSSP_DIR / f"{pid}_dssp.tsv"
        if not dssp.exists():
            fail += 1; fails.append((pid, "no_dssp")); continue
        sc = residue_scores(dssp)
        if not sc:
            fail += 1; fails.append((pid, "empty")); continue
        # cache par résidu
        cache = OUT_DIR / f"{pid}_a3d.csv"
        with cache.open("w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["residue_idx", "a3d_proxy"])
            for idx, s in enumerate(sc):
                w.writerow([idx, f"{s:.5f}"])
        rows.append((pid, mean(sc)))
        ok += 1

    with OUT_TSV.open("w") as fh:
        fh.write("pdb_chain\taggregation_score\n")
        for pid, s in rows:
            fh.write(f"{pid}\t{s:.5f}\n")

    print(f"[OK] {ok} réussis, {fail} échoués -> {OUT_TSV}")
    for pid, r in fails[:20]:
        print(f"  [FAIL] {pid} ({r})")
    if rows:
        vals = sorted(s for _, s in rows)
        n = len(vals)
        print(f"\nDistribution aggregation_score (n={n}) :")
        print(f"  min={vals[0]:.4f}  médiane={vals[n//2]:.4f}  max={vals[-1]:.4f}")
        print(f"  moyenne={mean(vals):.4f}")


if __name__ == "__main__":
    main()
