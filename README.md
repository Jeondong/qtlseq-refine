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
│   ├── qHT1_Chr1_p99.xlsx
│   ├── qHT1_Chr1_p99_sliding_window_multi.xlsx
│   ├── lesv_Chr11_p99.xlsx
│   └── lesv_Chr11_p99_sliding_window_multi.xlsx
│
└── README.md
```

---

## Data Files

All empirical example data files are located in the `data/` directory. Input files contain p99-filtered ΔSNP-index values extracted from primary QTL-seq analysis. Output files contain multi-scale sliding-window results generated from the corresponding input files.

| File                                        | Dataset              | Chromosome | Type                               |
| ------------------------------------------- | -------------------- | ---------- | ---------------------------------- |
| `abci8_Chr11_p99.xlsx`                      | OsABCI8              | Chr11      | Raw p99-filtered ΔSNP-index SNPs   |
| `abci8_Chr11_p99_sliding_window_multi.xlsx` | OsABCI8              | Chr11      | Multi-scale sliding-window results |
| `qHT1_Chr1_p99.xlsx`                        | qHT1 / pri-miR156b/c | Chr1       | Raw p99-filtered ΔSNP-index SNPs   |
| `qHT1_Chr1_p99_sliding_window_multi.xlsx`   | qHT1 / pri-miR156b/c | Chr1       | Multi-scale sliding-window results |
| `lesv_Chr11_p99.xlsx`                       | OsLESV               | Chr11      | Raw p99-filtered ΔSNP-index SNPs   |
| `lesv_Chr11_p99_sliding_window_multi.xlsx`  | OsLESV               | Chr11      | Multi-scale sliding-window results |
