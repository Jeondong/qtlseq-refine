# Supplementary Data S4

This directory contains all 11 binary parts of `Supplementary_Data_S4_snp_data.tsv.xz` (6,000,000 SNP records). Download the complete repository and run from its root:

```bash
python scripts/assemble_raw_data.py --extract
```

This joins the parts, verifies their SHA-256 hashes, verifies the reconstructed archive and decompressed TSV against the original submission data, and writes the archives and TSV files under `outputs/raw_supplementary/`. No external download or third-party Python package is needed. Individual part files are not standalone archives.

Only compression and packaging have changed. The decompressed TSV is byte-for-byte identical to the original `.tsv.gz` submission file. See [archive provenance](../../../provenance/raw_data_archives.json) for hashes and counts.
