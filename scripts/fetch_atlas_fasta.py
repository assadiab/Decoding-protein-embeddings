"""
fetch_atlas_fasta.py
Fetches FASTA sequences from the RCSB PDB REST API for all proteins in the
ATLAS id files (id_train.txt + id_test.txt).

Input IDs format: {pdb_code}_{chain}  (e.g. 6dnm_A)
Output: Datasets/ATLAS/atlas_sequences.fasta  (one entry per chain)

Usage:
    python scripts/fetch_atlas_fasta.py \
        --id-files deciphering/id_train.txt deciphering/id_test.txt \
        --out Datasets/ATLAS/atlas_sequences.fasta
"""

import argparse
import time
import requests
import os
from pathlib import Path


RCSB_FASTA_URL = "https://www.rcsb.org/fasta/entry/{pdb_id}"
DELAY_SEC = 0.2   # polite rate limiting: 5 req/s max


def load_ids(*id_files):
    ids = []
    for path in id_files:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line:
                    ids.append(line)
    return ids


def parse_chain_from_fasta(fasta_text, target_chain):
    """
    Parse RCSB FASTA response and return the sequence for the requested chain.
    RCSB header format:
      >4HHB_1|Chains A,B[auth A,B]|Hemoglobin subunit alpha|Homo sapiens (9606)
    We match chains by looking for 'auth {chain}' or 'Chains {chain}' patterns.
    """
    current_header = None
    current_seq_lines = []
    results = []

    for line in fasta_text.splitlines():
        if line.startswith(">"):
            if current_header is not None:
                results.append((current_header, "".join(current_seq_lines)))
            current_header = line
            current_seq_lines = []
        else:
            current_seq_lines.append(line.strip())
    if current_header is not None:
        results.append((current_header, "".join(current_seq_lines)))

    for header, seq in results:
        # Match chain: look for "[auth A]" or "[auth A,B]" or "Chains A"
        # RCSB format: "|Chains A[auth A]|" or "|Chains A,B[auth A,B]|"
        chain_upper = target_chain.upper()
        # Check [auth ...] section
        if f"[auth {chain_upper}]" in header or f"[auth {chain_upper}," in header \
                or f",{chain_upper}]" in header:
            return header, seq
        # Fallback: check "Chains X" without auth (some entries)
        import re
        chains_match = re.search(r"Chains? ([^|]+)", header)
        if chains_match:
            chains_str = chains_match.group(1)
            # Remove [auth ...] parts for simpler parsing
            chains_clean = re.sub(r"\[auth [^\]]+\]", "", chains_str)
            chain_list = [c.strip() for c in chains_clean.split(",")]
            if chain_upper in chain_list:
                return header, seq

    return None, None


def fetch_sequence(pdb_id, chain, retries=3):
    url = RCSB_FASTA_URL.format(pdb_id=pdb_id.upper())
    for attempt in range(retries):
        try:
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                header, seq = parse_chain_from_fasta(resp.text, chain)
                return header, seq
            elif resp.status_code == 404:
                print(f"  [404] {pdb_id} not found in RCSB")
                return None, None
            else:
                print(f"  [HTTP {resp.status_code}] {pdb_id} — retry {attempt+1}/{retries}")
                time.sleep(1)
        except requests.RequestException as e:
            print(f"  [ERR] {pdb_id}: {e} — retry {attempt+1}/{retries}")
            time.sleep(2)
    return None, None


def main():
    parser = argparse.ArgumentParser(description="Fetch ATLAS protein FASTA sequences from RCSB PDB")
    parser.add_argument("--id-files", nargs="+", required=True,
                        help="Paths to id_train.txt / id_test.txt")
    parser.add_argument("--out", required=True,
                        help="Output FASTA file path")
    parser.add_argument("--subset", type=int, default=None,
                        help="Only fetch the first N proteins (for testing)")
    args = parser.parse_args()

    ids = load_ids(*args.id_files)
    # Deduplicate while preserving order
    seen = set()
    unique_ids = [x for x in ids if not (x in seen or seen.add(x))]
    print(f"Total unique protein chains: {len(unique_ids)}")

    if args.subset:
        unique_ids = unique_ids[:args.subset]
        print(f"Running on subset of {len(unique_ids)} proteins")

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    ok, skipped = 0, 0
    with open(args.out, "w") as out_fasta:
        for i, protein_id in enumerate(unique_ids, 1):
            parts = protein_id.split("_")
            if len(parts) != 2:
                print(f"  [SKIP] Unexpected ID format: {protein_id}")
                skipped += 1
                continue
            pdb_code, chain = parts[0], parts[1]

            print(f"[{i:4d}/{len(unique_ids)}] {protein_id} ...", end=" ", flush=True)
            header, seq = fetch_sequence(pdb_code, chain)

            if seq:
                # Write with clean label matching PLM-API expectations
                out_fasta.write(f">{protein_id}\n{seq}\n")
                print(f"OK ({len(seq)} aa)")
                ok += 1
            else:
                print(f"SKIPPED (chain {chain} not found in entry {pdb_code.upper()})")
                skipped += 1

            time.sleep(DELAY_SEC)

    print(f"\nDone: {ok} fetched, {skipped} skipped → {args.out}")


if __name__ == "__main__":
    main()
