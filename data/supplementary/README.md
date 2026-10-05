# Supplementary data index

| Item | File or directory | What it contains |
| --- | --- | --- |
| Data S1 | [Supplementary_Data_S1_threshold_sim.tsv](Supplementary_Data_S1_threshold_sim.tsv) | 600 observed-percentile experiments; 5,400 result rows |
| Data S2 | [Supplementary_Data_S2_peakshift_sim.tsv](Supplementary_Data_S2_peakshift_sim.tsv) | 500 coordinate-assignment experiments; 2,500 result rows |
| Data S3 | `Supplementary_Data_S3_snp_data.tsv.gz` — generated as below | Raw simulated SNP draws for S1; original submission file about 61 MB |
| Data S4 | `Supplementary_Data_S4_snp_data.tsv.gz` — generated as below | Raw simulated SNP draws for S2; original submission file about 51 MB |
| Data S5 | [Supplementary_Data_S5/](Supplementary_Data_S5/) | F2 protocol, replicate records, thresholds and full/selected summaries |
| Data S6 | [Supplementary_Data_S6/](Supplementary_Data_S6/) plus [six original Excel files in the parent directory](../) | Rice support/extremum summaries and supported-window records |

S1/S2/S5 and the S6 summary were relocated from `results/` without changing their contents. S6 additionally includes `empirical_supported_windows.tsv`, a recalculation of all 15 empirical scans with empty windows omitted. Original six Excel files remain unchanged. Source hashes and verification records are in [../../provenance/](../../provenance/).

## Using Data S6

`empirical_supported_windows.tsv` contains 16,875 rows and columns `dataset`, `window`, `chrom`, `start`, `end`, `mid`, `value_mean`, and `count`. Coordinates are bp and windows are `[start, end)`. Every row has at least one SNP. Genuine zero means with positive count are retained. No significance test or upstream p99 selection is repeated by this export.

`empirical_summary.tsv` provides marker support and global extrema for each dataset/scale, including both coordinate conventions. The six original Excel files are included by reference to their stable parent-directory paths rather than duplicated here. See [the main README](../../README.md) for definitions and limitations.

## Generate large raw Data S3/S4

The original gzip files remain in the manuscript submission package and are not stored in this Git repository. Run these from the repository root to generate their simulated records:

```bash
python scripts/sim_threshold.py --raw-snps --output-dir outputs/descriptive
python scripts/peak_shift_sim.py --raw-snps --output-dir outputs/descriptive
```

Raw positions and ΔSNP-index are rounded to one and six decimal places. Gzip metadata and serialization can differ from the submission files. Default seeds and full-precision calculations reproduce S1/S2 metrics.

Tables S1–S5 are in [../../tables/](../../tables/). Table S5 describes empirical data; Data S5 is the F2 simulation bundle. Figures S1/S2 are marker support/population sensitivity, respectively.
