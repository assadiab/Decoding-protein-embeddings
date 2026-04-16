"""
download_atlas_data.py
Télécharge les données ATLAS (RMSF, Bfactor, Neq, PDB) depuis l'API ATLAS
pour toutes les protéines listées dans id_train.txt / id_test.txt.

Stratégie : télécharge le ZIP complet, extrait seulement les fichiers utiles
(TSV + PDB), puis supprime le ZIP. Résultat : ~200KB par protéine.

Usage :
    # Subset test (20 premières protéines) :
    python scripts/download_atlas_data.py --subset 20

    # Complet (1390 protéines, ~40 min en parallèle) :
    python scripts/download_atlas_data.py

    # Avec un seul fichier d'IDs :
    python scripts/download_atlas_data.py --id-files deciphering/id_train.txt
"""

import argparse
import io
import os
import time
import zipfile
import requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

ATLAS_API = "https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/analysis/{pdb_chain}"
# Fichiers à garder (les .xtc et .tpr font 40+ MB — on ne les veut pas)
KEEP_SUFFIXES = ["_RMSF.tsv", "_Bfactor.tsv", "_Neq.tsv", ".pdb", "_corresp.tsv"]
DEFAULT_ID_FILES = ["deciphering/id_train.txt", "deciphering/id_test.txt"]
DEFAULT_OUT_DIR = "Datasets/ATLAS/data"


def load_ids(*id_files):
    ids = []
    for path in id_files:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line:
                    ids.append(line)
    seen = set()
    return [x for x in ids if not (x in seen or seen.add(x))]


def download_protein(pdb_chain: str, out_dir: Path, max_retries: int = 3) -> tuple[str, bool, str]:
    """Download and extract relevant files for one protein. Returns (id, success, msg)."""
    # Check if already done (all expected files present)
    expected = [out_dir / f"{pdb_chain}{suf}" for suf in KEEP_SUFFIXES if suf != "_corresp.tsv"]
    if all(f.exists() for f in expected):
        return pdb_chain, True, "already exists"

    url = f"https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/analysis/{pdb_chain}"

    for attempt in range(max_retries):
        try:
            resp = requests.get(url, timeout=120, stream=True)
            if resp.status_code == 404:
                return pdb_chain, False, "404 Not Found"
            if resp.status_code != 200:
                if attempt < max_retries - 1:
                    time.sleep(2 ** attempt)
                    continue
                return pdb_chain, False, f"HTTP {resp.status_code}"

            # Load ZIP from memory
            content = resp.content
            with zipfile.ZipFile(io.BytesIO(content)) as zf:
                for name in zf.namelist():
                    for suffix in KEEP_SUFFIXES:
                        if name.endswith(suffix):
                            target = out_dir / name
                            target.write_bytes(zf.read(name))
                            break

            return pdb_chain, True, f"OK ({len(content) // 1024} KB downloaded)"

        except requests.RequestException as e:
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)
            else:
                return pdb_chain, False, f"Error: {e}"
        except zipfile.BadZipFile as e:
            return pdb_chain, False, f"Bad ZIP: {e}"

    return pdb_chain, False, "Max retries exceeded"


def main():
    parser = argparse.ArgumentParser(description="Download ATLAS MD data (RMSF, Bfactor, Neq, PDB)")
    parser.add_argument("--id-files", nargs="+", default=DEFAULT_ID_FILES,
                        help="Paths to id_train.txt / id_test.txt")
    parser.add_argument("--out", default=DEFAULT_OUT_DIR,
                        help="Output directory for extracted files")
    parser.add_argument("--subset", type=int, default=None,
                        help="Only download first N proteins (for testing)")
    parser.add_argument("--workers", type=int, default=4,
                        help="Number of parallel download threads (default: 4)")
    args = parser.parse_args()

    ids = load_ids(*args.id_files)
    if args.subset:
        ids = ids[:args.subset]

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Downloading ATLAS data for {len(ids)} proteins → {out_dir}")
    print(f"Parallel workers: {args.workers}")
    print(f"Files to extract: {KEEP_SUFFIXES}\n")

    ok, failed = 0, []

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(download_protein, pid, out_dir): pid for pid in ids}
        for i, future in enumerate(as_completed(futures), 1):
            pid, success, msg = future.result()
            status = "✓" if success else "✗"
            print(f"[{i:4d}/{len(ids)}] {status} {pid:12s}  {msg}")
            if success:
                ok += 1
            else:
                failed.append((pid, msg))

    print(f"\n=== Done: {ok}/{len(ids)} OK ===")
    if failed:
        print(f"Failed ({len(failed)}):")
        for pid, msg in failed:
            print(f"  {pid}: {msg}")
        fail_log = out_dir / "download_failures.txt"
        fail_log.write_text("\n".join(f"{pid}\t{msg}" for pid, msg in failed))
        print(f"Failures saved to {fail_log}")


if __name__ == "__main__":
    main()
