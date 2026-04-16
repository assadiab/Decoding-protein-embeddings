import os
import sys
import argparse
import torch
import numpy as np
import pandas as pd
from safetensors.torch import load_file

"""
Adapted version of full_prot_emb.py for 2nd-generation PLMs (ESM-C, Ankh2).

Key differences from the original:
  - Reads embeddings from a single .safetensors file (PLM-API output format)
    instead of per-protein .pt files.
    .safetensors layout: { "protein_id": tensor([seq_len, dim]), ... }
  - Supports new model types: esmc_300m (960), esmc_600m (1152), ankh2_large (1536)
  - Keeps all original 1G model types for baseline comparison.

Usage:
  python full_prot_emb_2g.py id_train.txt /path/to/esmc_300M.safetensors \\
      -emb esmc_300m -dssp_dir /path/to/dssp/ -md_dir /path/to/md/ --y acc
"""

# CONSTANTS

MAX_ACC = {
    'A': 106, 'C': 135, 'D': 163, 'E': 194, 'F': 197,
    'G': 84,  'H': 184, 'I': 169, 'K': 205, 'L': 164,
    'M': 188, 'N': 157, 'P': 136, 'Q': 198, 'R': 248,
    'S': 130, 'T': 142, 'V': 142, 'W': 227, 'Y': 222, '!': 1
}

DSSP_METRICS = ['acc', 'sec3', 'sec8', 'aa']
DYNAMICS_METRICS = ['rmsf', 'bfact', 'neq', 'bfact_norm']
MEDIANS_MD_DATA = {'rmsf': 0.106, 'bfact': 19.64, 'bfact_norm': -0.257, 'neq': 1}

# Embedding dimensions (1G + 2G)
EMB_DIMS = {
    # 1G baselines
    'esm6':   320,
    'esm12':  480,
    'esm30':  640,
    'esm33':  1280,
    'esm36':  2560,
    'esm48':  5120,
    'bert':   1024,
    'ankhb':  768,
    'ankhl':  1536,
    't5':     1024,
    'msa':    768,
    # 2G new models
    'esmc_300m':   960,
    'esmc_600m':   1152,
    'ankh2_large': 1536,
}

# Argument parsing
parser = argparse.ArgumentParser(
    description='Build per-residue dataset from PLM-API .safetensors embeddings + ATLAS labels.'
)
parser.add_argument('id_filename', type=str,
                    help='Path to id file (id_train.txt or id_test.txt)')
parser.add_argument('embed_file', type=str,
                    help='Path to the .safetensors embedding file (one file per model, all proteins)')
parser.add_argument('--y', choices=DSSP_METRICS + DYNAMICS_METRICS, default='acc',
                    help='Target variable Y')
parser.add_argument('-emb', choices=list(EMB_DIMS.keys()), default='esmc_300m',
                    help='Embedding model type')
parser.add_argument('-dssp_dir', type=str, required=True,
                    help='Directory containing {protein_id}_dssp.tsv files')
parser.add_argument('-md_dir', type=str, required=True,
                    help='Directory containing {protein_id}_prod_R1_merged.tsv files')
args = parser.parse_args()

print("ARGS", args)
print("Embedding type:", args.emb, f"({EMB_DIMS[args.emb]} dims)")

# Determine train/test split from filename
if "train" in args.id_filename:
    TRAIN_OR_TEST = "train"
elif "test" in args.id_filename:
    TRAIN_OR_TEST = "test"
else:
    print("Error: Neither 'train' nor 'test' in id filename. Exiting.")
    sys.exit(1)

# Load protein ID list
prot_ids = []
with open(args.id_filename, 'r') as f:
    for line in f:
        prot_ids.append(line.strip())
print(f"Loaded {len(prot_ids)} protein IDs from {args.id_filename}")

# Load all embeddings from .safetensors (lazy: only reads metadata until accessed)
print(f"Loading embeddings from {args.embed_file} ...")
all_embeddings = load_file(args.embed_file)
available_keys = set(all_embeddings.keys())
print(f"  → {len(available_keys)} proteins in safetensors file")

# Output directory
out_dir = f"./datasets_emb_{args.emb}"
os.makedirs(out_dir, exist_ok=True)
dataset_path = out_dir + "/"

# Array accumulators
concatenated_embeddings = None
Y_concat_continue = None
Y_concat_categ = None
prot_res_labels = None

skipped = 0
processed = 0

for prot_id in prot_ids:
    if prot_id not in available_keys:
        print(f"  [SKIP] {prot_id}: not found in safetensors file")
        skipped += 1
        continue

    # Load embedding tensor → numpy [seq_len, dim]
    emb = all_embeddings[prot_id]
    if isinstance(emb, torch.Tensor):
        emb = emb.numpy()

    # Load DSSP and MD annotation files
    dssp_path = os.path.join(args.dssp_dir, prot_id + '_dssp.tsv')
    md_path   = os.path.join(args.md_dir,   prot_id + '_prod_R1_merged.tsv')

    try:
        df_dssp = pd.read_csv(dssp_path, sep="\t") if args.y in DSSP_METRICS else None
        df_dyn  = pd.read_csv(md_path,   sep="\t") if args.y in DYNAMICS_METRICS else None
    except FileNotFoundError as e:
        print(f"  [SKIP] {prot_id}: {e}")
        skipped += 1
        continue

    # Validate sequence length matches embedding
    if args.y in DSSP_METRICS:
        if len(df_dssp['seq']) != emb.shape[0]:
            print(f"  [SKIP] {prot_id}: DSSP length {len(df_dssp['seq'])} != emb {emb.shape[0]}")
            skipped += 1
            continue
    else:  # DYNAMICS_METRICS
        if len(df_dyn['fasta_seq']) != emb.shape[0]:
            print(f"  [SKIP] {prot_id}: MD length {len(df_dyn['fasta_seq'])} != emb {emb.shape[0]}")
            skipped += 1
            continue

    print(f"  [OK] {prot_id}  shape={emb.shape}")
    processed += 1

    # Accumulate embeddings
    concatenated_embeddings = emb if concatenated_embeddings is None \
        else np.vstack((concatenated_embeddings, emb))

    # Build labels (prot_id_residueAA format)
    if args.y in DSSP_METRICS:
        labels = df_dssp['seq']
        res_labels = prot_id + "_" + labels
        prot_res_labels = res_labels if prot_res_labels is None \
            else np.hstack((prot_res_labels, res_labels))

        if args.y == 'acc':
            relative_sa = np.array(
                [(sa / MAX_ACC[aa]) * 100 for aa, sa in zip(labels, df_dssp['acc'])]
            )
            Y_concat_continue = relative_sa if Y_concat_continue is None \
                else np.hstack((Y_concat_continue, relative_sa))
            y_cat = np.where(relative_sa < 16, 0, 1)
            Y_concat_categ = y_cat if Y_concat_categ is None \
                else np.hstack((Y_concat_categ, y_cat))

        elif args.y == 'sec3':
            y_cat = df_dssp['DSSP_3'].values
            Y_concat_categ = y_cat if Y_concat_categ is None \
                else np.hstack((Y_concat_categ, y_cat))

        elif args.y == 'sec8':
            y_cat = df_dssp['DSSP_8'].values
            Y_concat_categ = y_cat if Y_concat_categ is None \
                else np.hstack((Y_concat_categ, y_cat))

        elif args.y == 'aa':
            y_cat = labels.str[-1].values
            Y_concat_categ = y_cat if Y_concat_categ is None \
                else np.hstack((Y_concat_categ, y_cat))

    else:  # DYNAMICS_METRICS
        labels = df_dyn['fasta_seq']
        res_labels = prot_id + "_" + labels
        prot_res_labels = res_labels if prot_res_labels is None \
            else np.hstack((prot_res_labels, res_labels))

        col_map = {'rmsf': 'RMSF', 'bfact': 'Bfactor', 'bfact_norm': 'BfactorZscore', 'neq': 'Neq'}
        y = df_dyn[col_map[args.y]].values
        Y_concat_continue = y if Y_concat_continue is None \
            else np.hstack((Y_concat_continue, y))

print(f"\nDone: {processed} proteins processed, {skipped} skipped.")
print(f"Embedding matrix shape: {concatenated_embeddings.shape}")

# Threshold continuous dynamics variables into binary categories
if args.y in DYNAMICS_METRICS:
    threshold = MEDIANS_MD_DATA[args.y]
    print(f"Binarizing {args.y} at median threshold = {threshold}")
    Y_concat_categ = np.where(Y_concat_continue <= threshold, 0, 1)

# Build column names
n_emb = EMB_DIMS[args.emb]
emb_cols = [f"emb{i}" for i in range(1, n_emb + 1)]

if Y_concat_categ is not None:
    concatenated_embeddings = np.concatenate(
        (concatenated_embeddings, Y_concat_categ.reshape(-1, 1)), axis=1
    )
if Y_concat_continue is not None:
    concatenated_embeddings = np.concatenate(
        (concatenated_embeddings, Y_concat_continue.reshape(-1, 1)), axis=1
    )

if args.y in ['aa', 'sec3', 'sec8']:
    columns = emb_cols + [args.y + '_categ']
else:
    columns = emb_cols + [args.y + '_categ', args.y + '_continue']

DF = pd.DataFrame(concatenated_embeddings, columns=columns, index=prot_res_labels)

out_file = dataset_path + f"emb_all_positions_{args.y}_{TRAIN_OR_TEST}.csv"
DF.to_csv(out_file)
print(f"\nDataset saved to: {out_file}")
print(f"Shape: {DF.shape}")
print(DF.head())
