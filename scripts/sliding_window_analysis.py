# Sliding Window Analysis Script
# Used for:
#  - lesv_Chr11_p99.xlsx
#  - abci8_Chr11_p99.xlsx
#  - abci8_Chr12_p99.xlsx
#
# Function:
#   Perform sliding window analysis on SNP position vs value (e.g., delta SNP-index)
#   and export results to a multi-sheet Excel file.
#
# Window settings:
#   1) 2 Mb / 100 Kb
#   2) 1 Mb / 50 Kb
#   3) 0.5 Mb / 10 Kb
#   4) 0.1 Mb / 5 Kb
#   5) 0.01 Mb / 5 Kb
#
# Output columns:
#   chrom, start, end, mid, value_mean, count

import pandas as pd

# ===============================
# USER SETTINGS
# ===============================
input_excel = "INPUT_FILE.xlsx"   # e.g. lesv_Chr11_p99.xlsx
chrom_name = "Chr11"              # change to Chr12 if needed
output_excel = "OUTPUT_sliding_window.xlsx"

# ===============================
# LOAD DATA
# ===============================
df = pd.read_excel(input_excel)

# Assume:
#   column 0 = genomic position (bp)
#   column 1 = value (e.g. delta SNP-index)
pos_col = df.columns[0]
val_col = df.columns[1]

# ===============================
# SLIDING WINDOW CONFIGS
# ===============================
configs = [
    ("2Mb_100Kb", int(2e6), int(100e3)),
    ("1Mb_50Kb",  int(1e6), int(50e3)),
    ("0.5Mb_10Kb",int(5e5), int(10e3)),
    ("0.1Mb_5Kb", int(1e5), int(5e3)),
    ("0.01Mb_5Kb",int(1e4), int(5e3)),
]

# ===============================
# RUN SLIDING WINDOW
# ===============================
max_pos = int(df[pos_col].max())
results = {}

for name, win, step in configs:
    rows = []
    for start in range(0, max_pos + 1, step):
        end = start + win
        sub = df[(df[pos_col] >= start) & (df[pos_col] < end)]
        if not sub.empty:
            rows.append({
                "chrom": chrom_name,
                "start": start,
                "end": end,
                "mid": (start + end) // 2,
                "value_mean": sub[val_col].mean(),
                "count": len(sub)
            })
    results[name] = pd.DataFrame(rows)

# ===============================
# SAVE TO EXCEL
# ===============================
with pd.ExcelWriter(output_excel, engine="xlsxwriter") as writer:
    for sheet, rdf in results.items():
        rdf.to_excel(writer, sheet_name=sheet, index=False)

print("Sliding window analysis completed.")
