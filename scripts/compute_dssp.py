"""
Compute DSSP for all PDB files in Datasets/ATLAS/data/ and write _dssp.tsv files.

DSSP v4 syntax: mkdssp input.pdb (no -i flag)
Handles GROMACS-style PDB headers by extracting ATOM/TER records first.

Output columns: seq, acc, DSSP_3, DSSP_8
  seq    : 1-letter amino acid
  acc    : solvent accessibility (integer, absolute Å²)
  DSSP_8 : 8-state (H G I E B T S or - for loop)
  DSSP_3 : 3-state (H=helix, E=strand, L=loop)
"""

import os
import subprocess
import tempfile
import pandas as pd
from pathlib import Path

ROOT = Path("/Volumes/T9_Assa/Cours/M2/S2/Projet Long/Code")
DATA_DIR = ROOT / "Datasets/ATLAS/data"
LIBCIFPP = str(ROOT / ".pixi/envs/default/share/libcifpp")

SS8_TO_SS3 = {
    'H': 'H', 'G': 'H', 'I': 'H',   # helices
    'E': 'E', 'B': 'E',               # strands
    'T': 'L', 'S': 'L', ' ': 'L', '-': 'L',  # loops
}


def clean_pdb(pdb_path: Path, tmp_path: str):
    """Extract ATOM/TER records — strips GROMACS headers and MODEL/ENDMDL."""
    with open(pdb_path) as f_in, open(tmp_path, 'w') as f_out:
        for line in f_in:
            if line.startswith(('ATOM', 'TER', 'END')):
                f_out.write(line)


def run_dssp(pdb_path: Path) -> str:
    """Run mkdssp on a PDB file (cleaned) and return classic DSSP output."""
    env = os.environ.copy()
    env['LIBCIFPP_DATA_DIR'] = LIBCIFPP

    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        tmp_path = tmp.name

    try:
        clean_pdb(pdb_path, tmp_path)
        result = subprocess.run(
            ['mkdssp', '--output-format', 'dssp', tmp_path],
            capture_output=True, text=True, env=env
        )
        return result.stdout
    finally:
        os.unlink(tmp_path)


def parse_dssp(dssp_output: str) -> pd.DataFrame:
    """Parse classic DSSP output into a DataFrame with seq/acc/DSSP_3/DSSP_8."""
    rows = []
    in_residues = False

    for line in dssp_output.splitlines():
        if line.startswith('  #  RESIDUE'):
            in_residues = True
            continue
        if not in_residues:
            continue

        # Skip chain breaks (marked with !)
        if len(line) < 38 or line[13] == '!':
            continue

        aa   = line[13]
        ss8  = line[16] if line[16] != ' ' else '-'
        ss3  = SS8_TO_SS3.get(line[16], 'L')
        acc  = int(line[35:38].strip())

        rows.append({'seq': aa, 'acc': acc, 'DSSP_8': ss8, 'DSSP_3': ss3})

    return pd.DataFrame(rows)


def main():
    pdb_files = sorted(DATA_DIR.glob('*.pdb'))
    print(f"Found {len(pdb_files)} PDB files")

    done = skipped = errors = 0

    for pdb_path in pdb_files:
        prot_id = pdb_path.stem
        out_path = DATA_DIR / f"{prot_id}_dssp.tsv"

        if out_path.exists() and out_path.stat().st_size > 0:
            skipped += 1
            continue

        try:
            dssp_out = run_dssp(pdb_path)
            df = parse_dssp(dssp_out)

            if df.empty:
                print(f"  [WARN] {prot_id}: empty DSSP output")
                errors += 1
                continue

            df.to_csv(out_path, sep='\t', index=False)
            done += 1

            if done % 100 == 0:
                print(f"  {done} done, {skipped} skipped, {errors} errors")

        except Exception as e:
            print(f"  [ERROR] {prot_id}: {e}")
            errors += 1

    print(f"\nDone: {done} generated, {skipped} already existed, {errors} errors")


if __name__ == '__main__':
    main()
