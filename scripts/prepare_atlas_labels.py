"""
prepare_atlas_labels.py
Fusionne les TSV bruts d'ATLAS (RMSF, Bfactor, Neq) en un fichier
`{prot_id}_prod_R1_merged.tsv` au format attendu par full_prot_emb_2g.py.

Format de sortie (colonnes) :
  fasta_seq | RMSF | Bfactor | BfactorZscore | Neq

- RMSF     = moyenne des 3 réplicats (RMSF_R1, RMSF_R2, RMSF_R3)
- Bfactor  = valeur expérimentale brute (NaN si absent → sera ignoré par le pipeline)
- BfactorZscore = (Bfactor - mean) / std, calculé sur les résidus avec B-factor valide
- Neq      = moyenne des 3 réplicats (Neq_R1, Neq_R2, Neq_R3)

Usage :
    # Traiter toutes les protéines dans Datasets/ATLAS/data/ :
    python scripts/prepare_atlas_labels.py

    # Tester sur un subset :
    python scripts/prepare_atlas_labels.py --subset 20
"""

import argparse
import os
import numpy as np
import pandas as pd
from pathlib import Path

DEFAULT_DATA_DIR = "Datasets/ATLAS/data"
DEFAULT_OUT_DIR = "Datasets/ATLAS/labels_md"


def prepare_one(prot_id: str, data_dir: Path, out_dir: Path) -> tuple[bool, str]:
    rmsf_path = data_dir / f"{prot_id}_RMSF.tsv"
    bfact_path = data_dir / f"{prot_id}_Bfactor.tsv"
    neq_path   = data_dir / f"{prot_id}_Neq.tsv"

    if not rmsf_path.exists():
        return False, f"RMSF file missing: {rmsf_path}"

    try:
        # RMSF: seq + 3 replicates → mean
        df_rmsf = pd.read_csv(rmsf_path, sep="\t")
        fasta_seq = df_rmsf["seq"].values
        rmsf_mean = df_rmsf[["RMSF_R1", "RMSF_R2", "RMSF_R3"]].mean(axis=1).values

        # Bfactor: single column, `""` → NaN
        bfact_vals = np.full(len(fasta_seq), np.nan)
        if bfact_path.exists():
            df_bfact = pd.read_csv(bfact_path, sep="\t", header=0,
                                   names=["Bfactor"])
            # Replace `""` string with NaN
            df_bfact["Bfactor"] = pd.to_numeric(df_bfact["Bfactor"], errors="coerce")
            if len(df_bfact) == len(fasta_seq):
                bfact_vals = df_bfact["Bfactor"].values
            else:
                # Length mismatch — align by position, fill rest with NaN
                n = min(len(df_bfact), len(fasta_seq))
                bfact_vals[:n] = df_bfact["Bfactor"].values[:n]

        # BfactorZscore = (B - mean) / std on valid residues only
        valid_mask = ~np.isnan(bfact_vals)
        bfact_zscore = np.full(len(fasta_seq), np.nan)
        if valid_mask.sum() > 1:
            mu = np.nanmean(bfact_vals)
            std = np.nanstd(bfact_vals)
            if std > 0:
                bfact_zscore[valid_mask] = (bfact_vals[valid_mask] - mu) / std

        # Neq: seq + 3 replicates → mean
        neq_mean = np.full(len(fasta_seq), np.nan)
        if neq_path.exists():
            df_neq = pd.read_csv(neq_path, sep="\t")
            if len(df_neq) == len(fasta_seq):
                neq_mean = df_neq[["Neq_R1", "Neq_R2", "Neq_R3"]].mean(axis=1).values

        # Build merged DataFrame
        df_out = pd.DataFrame({
            "fasta_seq":     fasta_seq,
            "RMSF":          rmsf_mean,
            "Bfactor":       bfact_vals,
            "BfactorZscore": bfact_zscore,
            "Neq":           neq_mean,
        })

        out_path = out_dir / f"{prot_id}_prod_R1_merged.tsv"
        df_out.to_csv(out_path, sep="\t", index=False)
        return True, f"{len(df_out)} residues"

    except Exception as e:
        return False, str(e)


def main():
    parser = argparse.ArgumentParser(
        description="Merge ATLAS TSV files into full_prot_emb_2g.py compatible format"
    )
    parser.add_argument("--data-dir", default=DEFAULT_DATA_DIR,
                        help="Directory containing downloaded ATLAS TSV files")
    parser.add_argument("--out", default=DEFAULT_OUT_DIR,
                        help="Output directory for merged TSV files")
    parser.add_argument("--subset", type=int, default=None,
                        help="Only process first N proteins")
    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    out_dir  = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Find all proteins in data_dir (look for _RMSF.tsv files)
    rmsf_files = sorted(data_dir.glob("*_RMSF.tsv"))
    prot_ids = [f.name.replace("_RMSF.tsv", "") for f in rmsf_files]

    if args.subset:
        prot_ids = prot_ids[:args.subset]

    print(f"Processing {len(prot_ids)} proteins → {out_dir}")

    ok, failed = 0, []
    for i, pid in enumerate(prot_ids, 1):
        success, msg = prepare_one(pid, data_dir, out_dir)
        status = "✓" if success else "✗"
        print(f"[{i:4d}/{len(prot_ids)}] {status} {pid:12s}  {msg}")
        if success:
            ok += 1
        else:
            failed.append((pid, msg))

    print(f"\n=== Done: {ok}/{len(prot_ids)} OK ===")
    if failed:
        print("Failed:", [pid for pid, _ in failed])


if __name__ == "__main__":
    main()
