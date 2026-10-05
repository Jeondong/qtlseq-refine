# Data

All Supplementary Data S1–S6, including the full raw simulated SNP records for S3/S4, are indexed in **[supplementary/](supplementary/)**. Tables remain in [../tables/](../tables/).

The six Excel files in this directory are preserved original empirical files and form part of Supplementary Data S6:

| Dataset | Retained SNP input | Original window output |
| --- | --- | --- |
| OsABCI8 | [abci8_Chr11_p99.xlsx](abci8_Chr11_p99.xlsx) | [abci8_Chr11_p99_sliding_window_multi.xlsx](abci8_Chr11_p99_sliding_window_multi.xlsx) |
| qHT1 | [qHT1_Chr1_p99.xlsx](qHT1_Chr1_p99.xlsx) | [qHT1_Chr1_p99_sliding_window_multi.xlsx](qHT1_Chr1_p99_sliding_window_multi.xlsx) |
| OsLESV | [lesv_Chr11_p99.xlsx](lesv_Chr11_p99.xlsx) | [lesv_Chr11_p99_sliding_window_multi.xlsx](lesv_Chr11_p99_sliding_window_multi.xlsx) |

The retained-input values and supported-window means/counts agree with recalculation. No numerical correction to those records was needed. Original outputs can contain empty windows encoded as mean zero: **exclude `count=0` rows before interpreting a profile.** Zero alone does not distinguish an empty window from a measured mean of zero.

For analysis, use the new [supported-window export](supplementary/Supplementary_Data_S6/empirical_supported_windows.tsv), which contains only nonempty windows, and the [verified empirical summary](supplementary/Supplementary_Data_S6/empirical_summary.tsv). These are derived from the original files, not newly collected biological data. Existing Excel paths remain stable.

From the repository root, `python scripts/verify_saved_results.py` verifies the original supported rows, both derived S6 files, and the F2 summaries. It writes fresh S6 exports under `outputs/verification/`.
