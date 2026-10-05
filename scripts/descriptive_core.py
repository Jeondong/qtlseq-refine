"""Gaussian signal-block model shared by the descriptive experiments.

Observed-percentile selection is not a null-calibrated significance test.
RandomState preserves the original NumPy MT19937 random stream.
"""
import numpy as np

CHROM_BP = 30_000_000
QTL_MB = 15.0
WINDOWS = {
    "2 Mb / 200 kb\n(conventional)": (2_000_000, 200_000),
    "500 kb / 50 kb\n(medium)": (500_000, 50_000),
    "100 kb / 10 kb\n(fine-scale)": (100_000, 10_000),
}

def make_chrom(rng):
    qpos = rng.uniform(13.5e6, 16.5e6, 1200)
    qval = np.clip(rng.normal(.72, .12, 1200), -1, 1)
    bpos = np.concatenate([rng.uniform(0, 13.5e6, 5400), rng.uniform(16.5e6, 30e6, 5400)])
    bval = np.clip(rng.normal(0, .20, 10800), -1, 1)
    pos = np.concatenate([qpos, bpos]); val = np.concatenate([qval, bval])
    order = np.argsort(pos)
    return pos[order], val[order]

def windows(pos, val, width, step, minimum=2, stop=CHROM_BP, include_empty=False):
    """Half-open windows; direct means preserve original summation/tie behavior."""
    starts = np.arange(0, stop, step)
    left = np.searchsorted(pos, starts); right = np.searchsorted(pos, starts + width)
    counts = right - left; valid = counts >= minimum
    means = np.full(len(starts), np.nan)
    pairs, inverse = np.unique(np.column_stack([left[valid], right[valid]]), axis=0, return_inverse=True)
    if len(pairs):
        means[valid] = np.array([val[a:b].mean() for a,b in pairs])[inverse]
    keep = np.ones(len(starts), dtype=bool) if include_empty else valid
    return starts[keep], starts[keep] + width / 2, means[keep], counts[keep]
