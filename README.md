# qtlseq-refine

**Multi-scale sliding-window refinement of QTL-seq significant ΔSNP-index signals**

This repository contains Python scripts and example data associated with the following manuscript:

> Jeon D, Shim K-C. *Refining QTL-seq Resolution by Re-analysis of Significant Delta SNP-Index Signals Using High-Resolution Sliding Windows.* (2026)

---

## Overview

Conventional QTL-seq analysis uses relatively large sliding windows, such as 2 Mb windows with 100 kb increments, to stabilize genome-wide SNP-index and ΔSNP-index estimates. Although this approach is robust for primary QTL detection, it often results in broad candidate intervals that contain many genes.

This repository provides scripts and example datasets for re-analyzing statistically significant ΔSNP-index signals using progressively smaller sliding-window parameters. The workflow is designed to refine QTL-seq signals within previously identified QTL regions without requiring additional sequencing or population development.

The repository includes:

1. A multi-scale sliding-window analysis script for empirical QTL-seq datasets
2. p99-filtered ΔSNP-index input data for three rice QTL-seq datasets:

   * OsABCI8
   * qHT1 / pri-miR156b/c
   * OsLESV
3. Monte Carlo simulation scripts for:

   * evaluating the p99 filtering threshold
   * verifying positional bias caused by sliding-window coordinate assignment

---

## Repository Structure

```text
qtlseq-refine/
│
├── scripts/
│   ├── sliding_window_analysis.py       # Empirical multi-scale sliding-window analysis
│   ├── sim_threshold.py                 # Monte Carlo simulation: p99 threshold evaluation
│   └── peak_shift_sim.py                # Monte Carlo simulation: positional bias verification
│
├── data/
│   ├── abci8_Chr11_p99.xlsx
│   ├── abci8_Chr11_p99_sliding_window_multi.xlsx
│   ├── abci8_Chr11_peak_5kb_2kb_sliding_window.xlsx
│   ├── qHT1_Chr1_p99.xlsx
│   ├── qHT1_Chr1_p99_sliding_window_multi.xlsx
│   ├── lesv_Chr11_p99.xlsx
│   └── lesv_Chr11_p99_sliding_window_multi.xlsx
│
└── README.md
```

---

## Scripts

### 1. `sliding_window_analysis.py`

Performs multi-scale sliding-window analysis on p99-filtered ΔSNP-index data.

**Input:**
Excel file containing at least the following two columns:

| Column           | Description                          |
| ---------------- | ------------------------------------ |
| `POSITION`       | Genomic position in base pairs       |
| `delta_SNPindex` | ΔSNP-index value after p99 filtering |

**Output:**
A multi-sheet Excel file containing sliding-window results for each window scale.

| Sheet name   | Window size  | Increment  |
| ------------ | ------------ | ---------- |
| `2Mb_100Kb`  | 2,000,000 bp | 100,000 bp |
| `1Mb_50Kb`   | 1,000,000 bp | 50,000 bp  |
| `0.5Mb_10Kb` | 500,000 bp   | 10,000 bp  |
| `0.1Mb_5Kb`  | 100,000 bp   | 5,000 bp   |
| `0.01Mb_5Kb` | 10,000 bp    | 5,000 bp   |

**Output columns:**

| Column       | Description                             |
| ------------ | --------------------------------------- |
| `chrom`      | Chromosome name                         |
| `start`      | Window start position in base pairs     |
| `end`        | Window end position in base pairs       |
| `mid`        | Window midpoint in base pairs           |
| `value_mean` | Mean ΔSNP-index value within the window |
| `count`      | Number of SNPs included in the window   |

**Usage:**

```bash
python scripts/sliding_window_analysis.py
```

Before running the script, edit the following variables inside `sliding_window_analysis.py` according to the target dataset:

```python
input_excel
chrom_name
output_excel
```

---

### 2. `sim_threshold.py`

Performs Monte Carlo simulations comparing three SNP filtering strategies:

* all SNPs
* p95-filtered SNPs
* p99-filtered SNPs

The simulation evaluates QTL signal enrichment, genome-wide noise, and signal discriminability across multiple sliding-window scales.

**Outputs:**

| Output file                         | Description                                       |
| ----------------------------------- | ------------------------------------------------- |
| `fig_main.png` / `fig_main.svg`     | 3 × 3 boxplot panel used for threshold comparison |
| `fig_traces.png` / `fig_traces.svg` | Representative 100 kb sliding-window profiles     |
| `table_threshold_clean.tsv`         | Summary statistics table                          |

**Key parameters:**

| Parameter        | Default | Description                                   |
| ---------------- | ------- | --------------------------------------------- |
| `CHROM_MB`       | 30      | Simulated chromosome length in Mb             |
| `QTL_POS_MB`     | 15.0    | True QTL position in Mb                       |
| `TOTAL_SNPS`     | 12,000  | Total number of SNPs per simulated chromosome |
| `QTL_MEAN`       | 0.72    | Mean ΔSNP-index value for QTL-linked SNPs     |
| `N_SIM`          | 600     | Number of Monte Carlo iterations              |
| `np.random.seed` | 42      | Random seed for reproducibility               |

**Usage:**

```bash
python scripts/sim_threshold.py
```

---

### 3. `peak_shift_sim.py`

Performs Monte Carlo simulations to evaluate positional bias caused by sliding-window coordinate assignment.

This script compares two coordinate-assignment methods:

1. **Start-position assignment**: the window value is plotted at the left boundary of the window.
2. **Center-position assignment**: the window value is plotted at the midpoint of the window.

The simulation tests whether the apparent downstream shift of ΔSNP-index peaks during progressive window-size reduction can be explained by coordinate assignment rather than biological displacement of the QTL signal.

**Outputs:**

| Output file                             | Description                                           |
| --------------------------------------- | ----------------------------------------------------- |
| `fig_N3_final.png` / `fig_N3_final.svg` | Six-panel figure summarizing positional bias analysis |
| Console summary table                   | Median positional bias by window scale                |

**Key parameters:**

| Parameter         | Default                      | Description                                   |
| ----------------- | ---------------------------- | --------------------------------------------- |
| `WINDOW_SIZES_MB` | `[2.0, 1.0, 0.5, 0.1, 0.01]` | Window sizes tested in Mb                     |
| `INCREMENT_RATIO` | 0.05                         | Increment size as a proportion of window size |
| `N_SIM`           | 500                          | Number of Monte Carlo iterations              |
| `np.random.seed`  | 2025                         | Random seed for reproducibility               |

**Usage:**

```bash
python scripts/peak_shift_sim.py
```

---

## Data Files

All empirical example data files are located in the `data/` directory. Input files contain p99-filtered ΔSNP-index values extracted from primary QTL-seq analysis. Output files contain multi-scale sliding-window results generated from the corresponding input files.

| File                                           | Dataset              | Chromosome | Type                                          |
| ---------------------------------------------- | -------------------- | ---------- | --------------------------------------------- |
| `abci8_Chr11_p99.xlsx`                         | OsABCI8              | Chr11      | Raw p99-filtered ΔSNP-index SNPs              |
| `abci8_Chr11_p99_sliding_window_multi.xlsx`    | OsABCI8              | Chr11      | Multi-scale sliding-window results            |
| `abci8_Chr11_peak_5kb_2kb_sliding_window.xlsx` | OsABCI8              | Chr11      | Fine-scale peak-region sliding-window results |
| `qHT1_Chr1_p99.xlsx`                           | qHT1 / pri-miR156b/c | Chr1       | Raw p99-filtered ΔSNP-index SNPs              |
| `qHT1_Chr1_p99_sliding_window_multi.xlsx`      | qHT1 / pri-miR156b/c | Chr1       | Multi-scale sliding-window results            |
| `lesv_Chr11_p99.xlsx`                          | OsLESV               | Chr11      | Raw p99-filtered ΔSNP-index SNPs              |
| `lesv_Chr11_p99_sliding_window_multi.xlsx`     | OsLESV               | Chr11      | Multi-scale sliding-window results            |

---

## Input File Format

The empirical input files named `*_p99.xlsx` contain p99-filtered ΔSNP-index values.

| Column           | Description                          |
| ---------------- | ------------------------------------ |
| `POSITION`       | Genomic position in base pairs       |
| `delta_SNPindex` | ΔSNP-index value after p99 filtering |

---

## Sliding-Window Output Format

The files named `*_sliding_window*.xlsx` contain sliding-window summary statistics.

| Column       | Description                             |
| ------------ | --------------------------------------- |
| `chrom`      | Chromosome name                         |
| `start`      | Window start position in base pairs     |
| `end`        | Window end position in base pairs       |
| `mid`        | Window midpoint in base pairs           |
| `value_mean` | Mean ΔSNP-index value within the window |
| `count`      | Number of SNPs included in the window   |

---

## Requirements

```text
python >= 3.8
numpy
pandas
matplotlib
openpyxl
xlsxwriter
```

Install dependencies using:

```bash
pip install numpy pandas matplotlib openpyxl xlsxwriter
```

---

## Reproducibility

The simulation scripts use fixed random seeds to ensure reproducibility.

| Script              | Random seed |
| ------------------- | ----------- |
| `sim_threshold.py`  | 42          |
| `peak_shift_sim.py` | 2025        |

Empirical sliding-window analysis can be reproduced by running `sliding_window_analysis.py` on each p99-filtered input file in the `data/` directory.

---

## Citation

If you use these scripts or data, please cite:

> Jeon D, Shim K-C. *Refining QTL-seq Resolution by Re-analysis of Significant Delta SNP-Index Signals Using High-Resolution Sliding Windows.* (2026)

---

## License

MIT License
