"""
QTL-seq threshold simulation - clean final version
논문 핵심 주장: p99 필터는 배경 SNP의 dilution effect를 제거하여
QTL 신호를 농축시킨다. 이는 작은 윈도우 사용을 가능하게 한다.

3가지 지표:
  A. Signal enrichment (QTL window mean Δ) - p99가 높을수록 좋음
  B. Noise level (background window std) - 전체 SNP 기반
  C. Discriminability (QTL signal / genome-wide std) - 실질 SNR
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

np.random.seed(42)

# ── Parameters (realistic bi-parental cross) ──────────────────────────────
CHROM_MB     = 30
QTL_POS_MB   = 15.0
QTL_HALF_MB  = 1.5
TOTAL_SNPS   = 12000        # ~400 SNPs/Mb
QTL_SNPS     = int(TOTAL_SNPS * 2 * QTL_HALF_MB / CHROM_MB)  # ~1200
BG_SNPS      = TOTAL_SNPS - QTL_SNPS

QTL_MEAN, QTL_STD = 0.72, 0.12
BG_STD            = 0.20
N_SIM             = 600

WINDOWS = {
    "2 Mb / 200 kb\n(conventional)": (2_000_000, 200_000),
    "500 kb / 50 kb\n(medium)":       (  500_000,  50_000),
    "100 kb / 10 kb\n(fine-scale)":   (  100_000,  10_000),
}
THRESHOLDS = ["All SNPs", "p95", "p99"]
COLORS = {"All SNPs": "#c0392b", "p95": "#e67e22", "p99": "#27ae60"}

def make_chrom():
    qtl_pos   = np.random.uniform((QTL_POS_MB-QTL_HALF_MB)*1e6,
                                   (QTL_POS_MB+QTL_HALF_MB)*1e6, QTL_SNPS)
    qtl_delta = np.clip(np.random.normal(QTL_MEAN, QTL_STD, QTL_SNPS), -1, 1)
    bg_left   = np.random.uniform(0, (QTL_POS_MB-QTL_HALF_MB)*1e6, BG_SNPS//2)
    bg_right  = np.random.uniform((QTL_POS_MB+QTL_HALF_MB)*1e6, CHROM_MB*1e6, BG_SNPS-BG_SNPS//2)
    bg_pos    = np.concatenate([bg_left, bg_right])
    bg_delta  = np.clip(np.random.normal(0, BG_STD, len(bg_pos)), -1, 1)
    pos   = np.concatenate([qtl_pos, bg_pos])
    delta = np.concatenate([qtl_delta, bg_delta])
    idx   = np.argsort(pos)
    return pos[idx], delta[idx]

def sliding_window(pos, val, win_bp, inc_bp):
    wp, wm, wn = [], [], []
    for s in np.arange(0, CHROM_MB*1e6, inc_bp):
        m = (pos >= s) & (pos < s+win_bp)
        if m.sum() >= 2:
            wp.append(s + win_bp/2)
            wm.append(val[m].mean())
            wn.append(m.sum())
    return np.array(wp), np.array(wm), np.array(wn)

records = []
example_traces = {wl: {t: [] for t in THRESHOLDS} for wl in WINDOWS}

for sim_i in range(N_SIM):
    pos, delta = make_chrom()
    cuts = {
        "All SNPs": 0.0,
        "p95":      np.percentile(np.abs(delta), 95),
        "p99":      np.percentile(np.abs(delta), 99),
    }

    for wl, (win_bp, inc_bp) in WINDOWS.items():
        for thr in THRESHOLDS:
            keep = np.abs(delta) >= cuts[thr]
            n_keep = keep.sum()
            if n_keep < 3: continue

            wp, wm, wn = sliding_window(pos[keep], delta[keep], win_bp, inc_bp)
            if len(wp) < 5: continue

            # Save example traces
            if sim_i < 5:
                example_traces[wl][thr].append((wp.copy(), wm.copy()))

            # 1. QTL signal: mean window value in QTL region
            qtl_m   = np.abs(wp/1e6 - QTL_POS_MB) <= 2.5
            qtl_sig = wm[qtl_m].mean() if qtl_m.sum() >= 2 else np.nan

            # 2. Genome-wide noise: std of ALL window means
            gw_noise = wm.std()

            # 3. Discriminability = qtl_sig / gw_noise
            disc = qtl_sig / gw_noise if (gw_noise > 0) else np.nan

            # 4. Peak detection
            peak_idx = np.argmax(wm)
            peak_pos = wp[peak_idx] / 1e6
            peak_err = abs(peak_pos - QTL_POS_MB)
            peak_val = wm[peak_idx]

            # 5. FWHM
            half  = peak_val / 2
            above = np.where(wm >= half)[0]
            fwhm  = (wp[above[-1]]-wp[above[0]])/1e6 if len(above)>=2 else np.nan

            records.append({
                "window": wl, "threshold": thr, "n_snps": n_keep,
                "qtl_signal": qtl_sig, "gw_noise": gw_noise,
                "discriminability": disc, "fwhm": fwhm,
                "peak_err": peak_err,
                "detected": int(peak_err <= 2.5),
            })

df = pd.DataFrame(records)

# ── Summary table ─────────────────────────────────────────────────────────
rows = []
for wl in WINDOWS:
    for thr in THRESHOLDS:
        s = df[(df["window"]==wl) & (df["threshold"]==thr)]
        rows.append({
            "Window":                       wl.replace("\n"," "),
            "Threshold":                    thr,
            "SNPs retained (mean)":         f"{s['n_snps'].mean():.0f}",
            "% of total":                   f"{s['n_snps'].mean()/TOTAL_SNPS*100:.1f}",
            "QTL signal Δ (median)":        f"{s['qtl_signal'].median():.3f}",
            "Genome-wide noise σ (median)": f"{s['gw_noise'].median():.4f}",
            "Discriminability (median)":    f"{s['discriminability'].median():.2f}",
            "Peak FWHM Mb (median)":        f"{s['fwhm'].median():.2f}",
            "Detection rate %":             f"{s['detected'].mean()*100:.1f}",
        })
df_table = pd.DataFrame(rows)
df_table.to_csv("/home/claude/table_threshold_clean.tsv", sep="\t", index=False)
print(df_table.to_string(index=False))

# ── Main figure ───────────────────────────────────────────────────────────
fig = plt.figure(figsize=(15, 11))
fig.patch.set_facecolor("#f7f8fa")
gs = GridSpec(3, 3, figure=fig, hspace=0.50, wspace=0.38)

win_list = list(WINDOWS.keys())

metric_info = [
    ("qtl_signal",       "QTL signal (mean ΔSNP-index)\nhigher = better signal enrichment", False, "ABC"),
    ("gw_noise",         "Genome-wide window noise (σ)\nlower = more stable estimation",    True,  "DEF"),
    ("discriminability", "Discriminability (signal / noise)\nhigher = clearer QTL peak",     False, "GHI"),
]

for row, (metric, ylabel, lower_better, panel_ids) in enumerate(metric_info):
    for col, wl in enumerate(win_list):
        ax = fig.add_subplot(gs[row, col])
        sub  = df[df["window"]==wl]
        data = [sub[sub["threshold"]==t][metric].dropna().values for t in THRESHOLDS]

        bp = ax.boxplot(data, patch_artist=True, widths=0.52,
                        medianprops=dict(color="white", linewidth=2.3),
                        flierprops=dict(marker=".", markersize=1.5, alpha=0.15),
                        whiskerprops=dict(color="#aaa", linewidth=0.9),
                        capprops=dict(color="#aaa", linewidth=0.9))
        for patch, thr in zip(bp["boxes"], THRESHOLDS):
            patch.set_facecolor(COLORS[thr]); patch.set_alpha(0.83)

        for i, (d, thr) in enumerate(zip(data, THRESHOLDS)):
            if len(d) > 0:
                med = np.median(d)
                yoff = abs(med) * 0.06 + 0.002
                va   = "bottom" if not lower_better else "top"
                sign = 1 if not lower_better else -1
                ax.text(i+1, med + sign*yoff, f"{med:.3f}" if metric=="gw_noise" else f"{med:.2f}",
                        ha="center", va=va, fontsize=7.5,
                        color=COLORS[thr], fontweight="bold")

        ax.set_xticks([1,2,3]); ax.set_xticklabels(THRESHOLDS, fontsize=9)
        win_short = wl.split("\n")[0]
        pl = panel_ids[col]
        ax.set_title(f"({pl})  {win_short}", fontsize=10, fontweight="bold")
        ax.set_ylabel(ylabel, fontsize=8.5)
        ax.grid(axis="y", alpha=0.22, linewidth=0.6)
        ax.set_facecolor("#fdfdfd")

        # Highlight best with coloured bottom spine
        medians = [np.median(d) if len(d)>0 else np.nan for d in data]
        valid = [(i,m) for i,m in enumerate(medians) if not np.isnan(m)]
        if valid:
            bi = min(valid,key=lambda x:x[1])[0] if lower_better \
                 else max(valid,key=lambda x:x[1])[0]
            ax.spines["bottom"].set_color(COLORS[THRESHOLDS[bi]])
            ax.spines["bottom"].set_linewidth(2.8)

legend_patches = [mpatches.Patch(facecolor=COLORS[t], label=t, alpha=0.85)
                  for t in THRESHOLDS]
fig.legend(handles=legend_patches, loc="lower center", ncol=3, fontsize=11,
           frameon=True, bbox_to_anchor=(0.5, -0.025),
           title="ΔSNP-index filtering threshold", title_fontsize=10)
fig.suptitle(
    "Simulation-based comparison of SNP filtering thresholds\n"
    "for high-resolution QTL-seq sliding-window re-analysis\n"
    f"(Monte Carlo: n = {N_SIM} iterations  |  chromosome = {CHROM_MB} Mb  |"
    f"  total SNPs = {TOTAL_SNPS:,})",
    fontsize=11, fontweight="bold", y=1.02)
plt.savefig("/home/claude/fig_main.png", dpi=180, bbox_inches="tight", facecolor="#f7f8fa")
plt.savefig("/home/claude/fig_main.svg", format="svg", bbox_inches="tight", facecolor="#f7f8fa")
plt.close()

# ── Example traces (fine-scale window) ────────────────────────────────────
fig2, axes2 = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
fig2.patch.set_facecolor("#f7f8fa")
fine_wl = "100 kb / 10 kb\n(fine-scale)"
for ax, thr in zip(axes2, THRESHOLDS):
    tr = example_traces[fine_wl][thr][:5]
    for wp, wm in tr:
        ax.plot(wp/1e6, wm, color=COLORS[thr], alpha=0.45, linewidth=1.0)
    # mean trace
    all_pos = np.arange(0, CHROM_MB, 0.1/10)  # dummy
    ax.axvline(QTL_POS_MB, color="#222", linewidth=1.4, linestyle="--", zorder=5)
    ax.axvspan(QTL_POS_MB-QTL_HALF_MB, QTL_POS_MB+QTL_HALF_MB,
               alpha=0.07, color="#222", zorder=1)
    n_median = int(df[(df["window"]==fine_wl)&(df["threshold"]==thr)]["n_snps"].median())
    ax.set_title(f"{thr}  (n ≈ {n_median} SNPs)\n100 kb window",
                 fontsize=10, fontweight="bold", color=COLORS[thr])
    ax.set_xlabel("Position (Mb)", fontsize=9)
    ax.set_xlim(0, CHROM_MB); ax.set_ylim(-0.15, 1.05)
    ax.grid(alpha=0.2); ax.set_facecolor("#fdfdfd")
axes2[0].set_ylabel("ΔSNP-index", fontsize=9)
fig2.suptitle(
    "Representative sliding-window profiles at fine scale (100 kb / 10 kb)\n"
    "Five simulation replicates per threshold condition  |  dashed line = true QTL position",
    fontsize=10, fontweight="bold")
plt.tight_layout()
plt.savefig("/home/claude/fig_traces.png", dpi=180, bbox_inches="tight", facecolor="#f7f8fa")
plt.savefig("/home/claude/fig_traces.svg", format="svg", bbox_inches="tight", facecolor="#f7f8fa")
plt.close()
print("\nDone.")
