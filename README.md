# qtlseq-refine

**Multi-scale sliding-window refinement of QTL-seq significant ΔSNP-index signals**

This repository contains Python scripts and example data associated with the following manuscript:

> Jeon D, Shim K-C. *Refining QTL-seq Resolution by Re-analysis of Significant Delta SNP-Index Signals Using High-Resolution Sliding Windows.* (2026)

---

## Overview

Conventional QTL-seq analysis uses large sliding windows (e.g., 2 Mb / 100 kb) to stabilize genome-wide SNP-index estimates, which limits mapping resolution. This repository provides:

1. A sliding window analysis script applied to empirical QTL-seq data (OsABCI8, OsLESV)
2. Simulation scripts validating (a) the p99 filtering threshold and (b) positional bias from window assignment convention

---

## Repository Structure

```
qtlseq-refine/
│
├── scripts/
│   ├── sliding_window_analysis.py      # Empirical sliding window (multi-scale)
│   ├── sim_threshold.py                # Monte Carlo: p99 threshold justification
│   └── peak_shift_sim.py               # Monte Carlo: positional bias verification
│
├── data/
│   ├── abci8_Chr11_p99.xlsx                        # OsABCI8 Chr11 p99-filtered SNPs
│   ├── abci8_Chr11_p99_sliding_window_multi.xlsx   # OsABCI8 Chr11 multi-scale SW results
│   ├── abci8_Chr11_peak_5kb_2kb_sliding_window.xlsx # OsABCI8 Chr11 fine-scale peak region
│   ├── lesv_Chr11_p99.xlsx                         # OsLESV Chr11 p99-filtered SNPs
│   └── lesv_Chr11_p99_sliding_window_multi.xlsx    # OsLESV Chr11 multi-scale SW results
│
└── README.md
```

---

## Scripts

### 1. `sliding_window_analysis.py`

Performs multi-scale sliding window analysis on p99-filtered ΔSNP-index data.

**Input:** Excel file with two columns — `POSITION` (bp) and `delta_SNPindex`

**Output:** Multi-sheet Excel file, one sheet per window scale

| Sheet name | Window size | Increment |
|------------|-------------|-----------|
| `2Mb_100Kb` | 2,000,000 bp | 100,000 bp |
| `1Mb_50Kb` | 1,000,000 bp | 50,000 bp |
| `0.5Mb_10Kb` | 500,000 bp | 10,000 bp |
| `0.1Mb_5Kb` | 100,000 bp | 5,000 bp |
| `0.01Mb_5Kb` | 10,000 bp | 5,000 bp |

**Output columns:** `chrom`, `start`, `end`, `mid`, `value_mean`, `count`

**Usage:**
```bash
python sliding_window_analysis.py
# Edit input_excel, chrom_name, output_excel inside the script before running
```

---

### 2. `sim_threshold.py`

Monte Carlo simulation (n = 600) comparing three SNP filtering strategies (All SNPs, p95, p99) across three window scales.

**Outputs:**
- `fig_main.png` / `fig_main.svg` — 3×3 boxplot panel (Figure 4 in manuscript)
- `fig_traces.png` / `fig_traces.svg` — Representative 100 kb sliding window profiles (Figure 3 in manuscript)
- `table_threshold_clean.tsv` — Summary statistics table

**Key parameters (editable in script):**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `CHROM_MB` | 30 | Simulated chromosome length (Mb) |
| `QTL_POS_MB` | 15.0 | True QTL position (Mb) |
| `TOTAL_SNPS` | 12,000 | Total SNPs on chromosome |
| `QTL_MEAN` | 0.72 | Mean ΔSNP-index at QTL locus |
| `N_SIM` | 600 | Number of Monte Carlo iterations |
| `np.random.seed` | 42 | Random seed (for reproducibility) |

**Usage:**
```bash
python sim_threshold.py
```

---

### 3. `peak_shift_sim.py`

Monte Carlo simulation (n = 500) verifying that the downstream peak shift observed with decreasing window size reflects window start vs. center position assignment convention, not biological signal displacement.

**Outputs:**
- `fig_N3_final.png` / `fig_N3_final.svg` — 6-panel figure (Figure 8 in manuscript)
- Console summary table of median positional bias per window scale

**Key parameters (editable in script):**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `WINDOW_SIZES_MB` | [2.0, 1.0, 0.5, 0.1, 0.01] | Window sizes tested (Mb) |
| `INCREMENT_RATIO` | 0.05 | Increment = window × 5% |
| `N_SIM` | 500 | Number of Monte Carlo iterations |
| `np.random.seed` | 2025 | Random seed (for reproducibility) |

**Usage:**
```bash
python peak_shift_sim.py
```

---

## Data Files

All data files are located in `data/`. Input files contain p99-filtered ΔSNP-index values extracted from primary QTL-seq analysis.

| File | Dataset | Chromosome | Type |
|------|---------|------------|------|
| `abci8_Chr11_p99.xlsx` | OsABCI8 | Chr11 | Raw p99-filtered SNPs |
| `abci8_Chr11_p99_sliding_window_multi.xlsx` | OsABCI8 | Chr11 | Multi-scale SW results |
| `abci8_Chr11_peak_5kb_2kb_sliding_window.xlsx` | OsABCI8 | Chr11 | Fine-scale peak region SW |
| `lesv_Chr11_p99.xlsx` | OsLESV | Chr11 | Raw p99-filtered SNPs |
| `lesv_Chr11_p99_sliding_window_multi.xlsx` | OsLESV | Chr11 | Multi-scale SW results |

### Input file format (`*_p99.xlsx`)

| Column | Description |
|--------|-------------|
| `POSITION` | Genomic position (bp) |
| `delta_SNPindex` | ΔSNP-index value (p99-filtered) |

### Sliding window output format (`*_sliding_window*.xlsx`)

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome name |
| `start` | Window start position (bp) |
| `end` | Window end position (bp) |
| `mid` | Window midpoint (bp) — used as x-axis for plotting |
| `value_mean` | Mean ΔSNP-index within window |
| `count` | Number of SNPs in window |

---

## Requirements

```
python >= 3.8
numpy
pandas
matplotlib
openpyxl
xlsxwriter
```

Install dependencies:
```bash
pip install numpy pandas matplotlib openpyxl xlsxwriter
```

---

## Citation

If you use these scripts or data, please cite:

> Jeon D, Shim K-C. *Refining QTL-seq Resolution by Re-analysis of Significant Delta SNP-Index Signals Using High-Resolution Sliding Windows.* (2026)

---

## License

MIT License
