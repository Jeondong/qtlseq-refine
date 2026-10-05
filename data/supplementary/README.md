# Supplementary data index

[Download the complete repository ZIP](https://github.com/Jeondong/qtlseq-refine/archive/refs/heads/main.zip) to obtain all data files and reconstruction scripts.

| Item | File or directory | What it contains |
| --- | --- | --- |
| Data S1 | [Supplementary_Data_S1_threshold_sim.tsv](Supplementary_Data_S1_threshold_sim.tsv) | 600 observed-percentile experiments; 5,400 result rows |
| Data S2 | [Supplementary_Data_S2_peakshift_sim.tsv](Supplementary_Data_S2_peakshift_sim.tsv) | 500 coordinate-assignment experiments; 2,500 result rows |
| Data S3 | [Supplementary_Data_S3/](Supplementary_Data_S3/) | Raw simulated SNP draws for S1; 7,200,000 SNP records; XZ archive 47.0 MB |
| Data S4 | [Supplementary_Data_S4/](Supplementary_Data_S4/) | Raw simulated SNP draws for S2; 6,000,000 SNP records; XZ archive 39.2 MB |
| Data S5 | [Supplementary_Data_S5/](Supplementary_Data_S5/) | F2 protocol, replicate records, thresholds and full/selected summaries |
| Data S6 | [Supplementary_Data_S6/](Supplementary_Data_S6/) plus [six original Excel files in the parent directory](../) | Rice support/extremum summaries and supported-window records |

S1/S2/S5 and the S6 summary were relocated from `results/` without changing their contents. S6 additionally includes `empirical_supported_windows.tsv`, a recalculation of all 15 empirical scans with empty windows omitted. Original six Excel files remain unchanged. Source hashes and verification records are in [../../provenance/](../../provenance/).

## Using Data S6

`empirical_supported_windows.tsv` contains 16,875 rows and columns `dataset`, `window`, `chrom`, `start`, `end`, `mid`, `value_mean`, and `count`. Coordinates are bp and windows are `[start, end)`. Every row has at least one SNP. Genuine zero means with positive count are retained. No significance test or upstream p99 selection is repeated by this export.

`empirical_summary.tsv` provides marker support and global extrema for each dataset/scale, including both coordinate conventions. The six original Excel files are included by reference to their stable parent-directory paths rather than duplicated here. See [the main README](../../README.md) for definitions and limitations.

## Download, verify and regenerate Data S3/S4

Both raw datasets are fully included above as binary parts of `.tsv.xz` archives (14 parts for S3; 11 for S4). Individual parts cannot be opened independently. After downloading the complete repository, run `python scripts/assemble_raw_data.py --extract` from its root to restore both complete archives and TSV files under `outputs/raw_supplementary/`. No external data download or third-party package is required. Only compression has changed: the decompressed TSV bytes match the original submission gzip files exactly. Archive SHA-256, decompressed TSV SHA-256, original gzip SHA-256 and record counts are in [raw_data_archives.json](../../provenance/raw_data_archives.json).

From the repository root, run `python scripts/verify_saved_results.py --verify-raw` to check both archives and their uncompressed contents. After assembly, Pandas reads the complete XZ archives directly with `pd.read_csv(path, sep="\t")`; Python's built-in `lzma` module can also decompress them.

To regenerate the simulated records as gzip files, run:

```bash
python scripts/sim_threshold.py --raw-snps --output-dir outputs/descriptive
python scripts/peak_shift_sim.py --raw-snps --output-dir outputs/descriptive
```

Raw positions and ΔSNP-index are rounded to one and six decimal places. Archive metadata and serialization can differ when regenerating the records. Default seeds and full-precision calculations reproduce S1/S2 metrics.

Tables S1–S5 are in [../../tables/](../../tables/). Table S5 describes empirical data; Data S5 is the F2 simulation bundle. Figures S1/S2 are marker support/population sensitivity, respectively.
