# qtlseq-refine

Multi-scale sliding-window summaries of retained ΔSNP-index signals, with descriptive simulations, an F2 population benchmark, and three rice examples.

Associated manuscript: Jeon D, Shim K-C. *Refining QTL-seq Resolution by Re-analysis of Significant Delta SNP-Index Signals Using High-Resolution Sliding Windows.* (2026).

## Revision status

The 2026-10-05 update incorporates the September simulation revision and verification of the supplied empirical files. The original six rice Excel files are unchanged. Existing F2 replicate results were reaggregated; no new F2 populations were generated for the reported revision. Previous repository versions remain available in Git history.

Filtering can narrow a plotted or called region while reducing true-position inclusion. A narrower region alone does not establish better localization. The rice examples illustrate retained-marker support and coordinate conventions; they do not independently establish mapping accuracy.

## Install and verify

Use Python 3.10 or later and run commands from the repository root:

```bash
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
python scripts/verify_saved_results.py
```

The verifier reaggregates 364,000 F2 replicate-method records (26,000 evaluation populations), checks 1,274 summary values, and recalculates 16,875 supported empirical windows across 15 settings. It writes a report under `outputs/verification/`. The additional 26,000 null calibration populations are documented by the saved protocol; their individual draws were not retained. See `provenance/` for file hashes and validation results. Add `--verify-raw` to also decompress and verify Data S3/S4 against the published archive and TSV hashes.

## Files and supplementary numbering

Start at the [data index](data/README.md) or [supplementary data index](data/supplementary/README.md). Saved supplementary results are now grouped under `data/supplementary/`; the eight files moved from `results/` retain their original contents.

| Location | Contents | Supplement |
| --- | --- | --- |
| `data/supplementary/Supplementary_Data_S1_threshold_sim.tsv` | 600 observed-percentile experiments, seed 42; 5,400 rows | Data S1 |
| `data/supplementary/Supplementary_Data_S2_peakshift_sim.tsv` | 500 coordinate experiments, seed 2025; 2,500 rows | Data S2 |
| [Data S3](data/supplementary/Supplementary_Data_S3/) and [Data S4](data/supplementary/Supplementary_Data_S4/) | Raw simulated SNP draws corresponding to S1 and S2, in split XZ archives | Data S3 and S4 |
| `data/supplementary/Supplementary_Data_S5/` | F2 protocol, thresholds, replicate records, full and selected summaries | Data S5 |
| `data/` and `data/supplementary/Supplementary_Data_S6/` | Six original rice Excel files and verified summaries | Data S6 |
| `tables/Supplementary_Tables.xlsx` and `tables/*.tsv` | Tables S1, S2, S3, S4, S5A and S5B | Tables S1–S5 |
| `figures/` | Revised PNG/SVG exports for Figures 1, 3, 4, 8, S1 and S2 | See below |

**Figure S1 is marker support; Figure S2 is population sensitivity.** The former `Figure_9_population_sensitivity` is now Figure S2. Publication exports retain their author-edited layout; the plotting script reproduces the reported quantities with potentially different typography/layout. Figures 5–7 are not regenerated here.

All six supplementary datasets are available through the [supplementary data index](data/supplementary/README.md). Data S3/S4 are published as lossless split XZ archives; all parts are included in the repository and `python scripts/assemble_raw_data.py --extract` restores both TSV files under `outputs/raw_supplementary/`. Their decompressed TSV bytes are identical to those in the original submission gzip files. They contain simulated SNPs, not raw sequencing reads or unfiltered empirical SNPs. “Table S5” and “Data S5” denote different supplementary items.

## Descriptive simulations and figures

```bash
python scripts/sim_threshold.py --output-dir outputs/descriptive
python scripts/peak_shift_sim.py --output-dir outputs/descriptive
python scripts/plot_revision_figures.py --output-dir outputs/figures
```

To regenerate raw simulated SNP draws (the scripts export gzip; the published archives use XZ):

```bash
python scripts/sim_threshold.py --raw-snps --output-dir outputs/descriptive
python scripts/peak_shift_sim.py --raw-snps --output-dir outputs/descriptive
```

Default seeds and replicate counts reproduce saved Data S1/S2 metrics. Archive metadata and floating-point text formatting may differ when rerunning the generators. Raw exports round positions to one decimal place and ΔSNP-index to six decimal places, matching the submission format; use the seeded generator for full-precision metric reproduction.

These experiments draw 12,000 SNPs over a 30 Mb chromosome with a signal centered at 15 Mb. Observed p95/p99 are percentiles of absolute values in each simulated chromosome, **not null-model significance thresholds**. Windows require two retained SNPs. Empirical and F2 support rules differ, as specified below.

Data S1 retains legacy column names:

- `qtl_signal`: mean of retained window means with midpoints within 2.5 Mb of 15 Mb.
- `gw_noise`: population SD of all retained window means, including signal-bearing windows.
- `discriminability`: signal-region mean divided by that SD; not an independent background-noise or localization metric.
- `peak_err`: midpoint error in Mb. `fwhm`: span between the outermost half-height window midpoints, including gaps; not a confidence interval.
- `detected`: peak error at most 2.5 Mb in this descriptive model.

In Data S2, midpoint assignment adds exactly half a window width to the start coordinate. This changes directional bias but leaves the spread of signed errors unchanged. Reduction in absolute median signed error and reduction in median absolute error are distinct statistics; Table S2 reports both. Figure 8 uses the first seed-2025 example, including its downstream 100 kb start-coordinate peak.

## Empirical rice scans

```bash
python scripts/sliding_window_analysis.py --input data/abci8_Chr11_p99.xlsx --chrom Chr11 --output-dir outputs/OsABCI8
python scripts/sliding_window_analysis.py --input data/qHT1_Chr1_p99.xlsx --chrom Chr1 --output-dir outputs/qHT1
python scripts/sliding_window_analysis.py --input data/lesv_Chr11_p99.xlsx --chrom Chr11 --output-dir outputs/OsLESV
```

The first two columns supply position in bp and retained ΔSNP-index. The same supplied input is scanned independently at each setting:

| Window | Step | Output stem |
| --- | --- | --- |
| 2 Mb | 100 kb | `2Mb_100Kb` |
| 1 Mb | 50 kb | `1Mb_50Kb` |
| 500 kb | 10 kb | `0.5Mb_10Kb` |
| 100 kb | 5 kb | `0.1Mb_5Kb` |
| 10 kb | 5 kb | `0.01Mb_5Kb` |

Each scan TSV contains `chrom`, `start`, `end`, `mid`, `value_mean`, and `count`. Windows are half-open `[start, end)` and `mid = start + width/2`. The default requires one retained SNP (`--min-snps 1`); unsupported windows are omitted. Original Excel outputs sometimes encode empty windows as mean zero; verification excludes `count=0` rows. The combined [supported-window export in Data S6](data/supplementary/Supplementary_Data_S6/empirical_supported_windows.tsv) adds `dataset` and `window` identifiers and preserves genuine zero means with positive SNP count. The verifier regenerates and checks all 16,875 rows.

Filenames preserve the original `p99` labels. The retained inputs do not include all upstream data/code needed to independently reconstruct selection. This scan is not an automated iterative interval-selection algorithm.

Inputs contain 24,443 SNPs for OsABCI8, 12,068 for qHT1, and 103 for OsLESV. At 10 kb, median support is 14, 24 and 1 SNPs; single-SNP windows comprise 12.45%, 4.84% and 97.5% of supported windows. Small windows can expose sparse support rather than add mapping information.

Table S5B selects the global minimum over supplied OsABCI8/qHT1 inputs and global maximum over OsLESV, taking the earliest start in an exact tie. Extrema are selected before referring to gene coordinates and can differ from local arrows in Figures 5–7. Distance is to the nearest locus boundary (zero inside); both start and midpoint distances are reported. The original empirical plotting convention still requires author confirmation.

## F2 population benchmark

Seed 20260914 is used for 13 controlled conditions, each with 2,000 null calibration, 1,000 independent null-test and 1,000 QTL-test populations (52,000 overall). The single-QTL F2 model uses 400 individuals and varies bulk size, depth, marker density, recombination map, effect size and error. These are sensitivity scenarios, not estimates fitted to rice data. Settings and limitations are in `data/supplementary/Supplementary_Data_S5/benchmark_protocol.json`.

Here p95/p99 are bulk-aware **null** cutoffs for individual absolute ΔSNP-index values. Window detection uses a separate held-out chromosome-maximum null threshold at nominal alpha 0.05, with strict exceedance. This concerns one chromosome, not genome-wide error. ALL/p95/p99 require ten retained SNPs; `p99_n1` is a one-SNP sensitivity analysis. Suffixes 2000/500/100 kb denote window width. The marker-fraction comparator uses the same p99 cutoff and is not PyBSASeq. The full summary preserves a reimplemented G-statistic comparator; no official comparator package execution is claimed.

Record and summary definitions:

- `stage` is `null` or `qtl`; `replicate` is paired across methods within scenario/stage.
- `detected`: a finite window strictly exceeds its threshold. `covered`: the true coordinate falls in the union of exceeding windows.
- `total_width_mb`: union length of complete exceeding windows; zero for nondetections. This is not a location confidence interval.
- `peak_error_mb`: midpoint error of the maximum finite score, including nondetections; ties within 1e-10 select the earliest window. Missing values mean no finite peak.
- Inclusion and no-finite-peak rates use all QTL tests; null detection uses all null tests. Median called width is conditional on detection; median peak error is conditional on a finite peak. Medians may use different populations.

At 100 kb, the ten-SNP p99 calibration threshold is negative infinity in 12/13 conditions because null profiles usually lack supported windows. Any finite test profile then exceeds that cutoff. This and differing realized null detection rates matter for comparisons.

At baseline 100 kb, ALL → p99 reduces median called width from 16.080 to 8.730 Mb, while median peak error stays 0.295 Mb and inclusion changes from 100.0% to 98.4%. Under weak effects, inclusion changes from 77.2% to 26.0%. This supports a conditional width/inclusion trade-off, not a universal accuracy gain.

Generate a complete new benchmark with:

```bash
python scripts/population_benchmark.py --n 1000 --n-cal 2000 --scenario all --output-dir outputs/f2_new --skip-plots
```

This is substantially more expensive than verification. Choose a new output directory for each run. Numerical model functions are retained from the source benchmark; this revision makes paths and command-line use portable. Saved records remain the source for manuscript tables and supplement figures.

## Citation and license

Please cite the manuscript above. The repository's stated license is MIT License.
