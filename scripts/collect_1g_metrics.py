"""Parse run_analyses_1g.sh log and emit results/dt_results_1g.csv."""
import re
import sys
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

EMB_TO_MODEL = {
    "esm33": "esm2_650m",
    "ankhl": "ankh_large",
    "t5": "prot_t5",
}

# Find most recent 1G analysis log
logs = sorted(RESULTS.glob("analyses_1g_*.log"))
if not logs:
    print("No analyses_1g_*.log found in results/")
    sys.exit(1)

log_path = logs[-1]
print(f"Parsing {log_path.name}")
text = log_path.read_text()

rows = []
header_re = re.compile(r"\[.*?\] === (\w+) × (\w+) \(dt\) ===")
f1_re = re.compile(r"F1 score:\s*([\d.]+)")
mcc_re = re.compile(r"MCC:([\d.\-]+)")

lines = text.splitlines()
i = 0
while i < len(lines):
    m = header_re.search(lines[i])
    if m:
        emb, task = m.group(1), m.group(2)
        model = EMB_TO_MODEL.get(emb, emb)
        f1_vals, mcc_vals = [], []
        for j in range(i + 1, min(i + 120, len(lines))):
            mf = f1_re.search(lines[j])
            mm = mcc_re.search(lines[j])
            if mf:
                f1_vals.append(float(mf.group(1)))
            if mm:
                mcc_vals.append(float(mm.group(1)))
            if "done" in lines[j] and emb in lines[j] and task in lines[j]:
                break
        if len(f1_vals) >= 2 and len(mcc_vals) >= 2:
            rows.append({"model": model, "task": task, "split": "train",
                         "F1": f1_vals[0], "MCC": mcc_vals[0]})
            rows.append({"model": model, "task": task, "split": "test",
                         "F1": f1_vals[1], "MCC": mcc_vals[1]})
    i += 1

if not rows:
    print("No metrics found — log may be incomplete.")
    sys.exit(1)

df = pd.DataFrame(rows)
out = RESULTS / "dt_results_1g.csv"
df.to_csv(out, index=False)
print(f"Wrote {len(rows)} rows to {out}")
print(df.pivot_table(index="task", columns="model", values="F1").round(3))
