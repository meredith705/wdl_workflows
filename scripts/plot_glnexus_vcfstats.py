#!/usr/bin/env python3
"""
plot_vcf_gq.py — Plot GQ distribution and per-sample variant counts
               from a GLnexus-merged DeepVariant VCF.

Usage:
    python plot_vcf_gq.py input.vcf[.gz] [-o output_prefix]
                          [--min-gq 0] [--pass-only] [--max-samples 200]

Outputs:
    <prefix>_gq_report.png  — two-panel figure (GQ histogram + per-sample violin)

Requirements:
    pip install cyvcf2 matplotlib seaborn pandas numpy
"""

import argparse
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

try:
    from cyvcf2 import VCF
except ImportError:
    sys.exit("ERROR: cyvcf2 not installed. Run: pip install cyvcf2")


# ── parse ─────────────────────────────────────────────────────────────────────

def parse_vcf(vcf_path: str, pass_only: bool, min_gq: float):
    """
    Returns:
        gq_all       : flat list of all non-missing GQ values across samples
        sample_counts: dict {sample_name: variant_count}  (GQ >= min_gq, GT called)
        sample_gqs   : dict {sample_name: [gq, ...]}
        n_variants   : total variant records parsed
    """
    vcf = VCF(vcf_path)
    samples = vcf.samples
    n_samples = len(samples)

    if n_samples == 0:
        sys.exit("ERROR: No samples found in VCF.")

    print(f"[…] {n_samples} samples detected")

    gq_all = []
    sample_counts = defaultdict(int)
    sample_gqs    = defaultdict(list)
    n_variants = 0

    for v in vcf:
        filt = v.FILTER
        filter_str = filt if filt else "PASS"
        if pass_only and filter_str != "PASS":
            continue

        n_variants += 1

        # GQ is a FORMAT field — one value per sample
        try:
            gqs = v.format("GQ")   # numpy array shape (n_samples, 1) or (n_samples,)
        except KeyError:
            continue

        gqs = gqs.flatten().astype(float)
        # GT array to detect missing / ref-only calls
        gts = v.gt_types  # 0=HOM_REF, 1=HET, 2=UNKNOWN, 3=HOM_ALT

        for i, (gq, gt) in enumerate(zip(gqs, gts)):
            if gq < 0 or np.isnan(gq):   # cyvcf2 encodes missing as -1
                continue
            if gt == 2:                   # unknown / missing GT
                continue
            if gq < min_gq:
                continue
            gq_all.append(gq)
            sample_gqs[samples[i]].append(gq)
            # Count as a variant call if not HOM_REF
            if gt != 0:
                sample_counts[samples[i]] += 1

    vcf.close()

    # Ensure every sample appears even if it had zero calls
    for s in samples:
        if s not in sample_counts:
            sample_counts[s] = 0

    print(f"[…] {n_variants:,} records parsed, {len(gq_all):,} genotype GQ values collected")
    return gq_all, sample_counts, sample_gqs, n_variants


# ── plot ──────────────────────────────────────────────────────────────────────

BG       = "#1a1a2e"
PANEL_BG = "#16213e"
ACCENT   = "#ffdd57"
TEXT     = "#eeeeff"
SPINE    = "#444466"

def style_ax(ax):
    ax.set_facecolor(PANEL_BG)
    for spine in ax.spines.values():
        spine.set_color(SPINE)
    ax.tick_params(colors=TEXT, labelsize=9)
    ax.xaxis.label.set_color(TEXT)
    ax.yaxis.label.set_color(TEXT)
    ax.title.set_color(TEXT)


def plot_report(gq_all, sample_counts, sample_gqs, n_variants,
                out_path, vcf_path, max_samples):

    n_samples = len(sample_counts)
    counts_series = pd.Series(sample_counts)

    fig = plt.figure(figsize=(18, 12))
    fig.patch.set_facecolor(BG)
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.32,
                           height_ratios=[1, 1])

    ax_gq_hist   = fig.add_subplot(gs[0, 0])
    ax_gq_violin = fig.add_subplot(gs[0, 1])
    ax_cnt_violin= fig.add_subplot(gs[1, :])   # full-width bottom panel

    for ax in [ax_gq_hist, ax_gq_violin, ax_cnt_violin]:
        style_ax(ax)

    # ── 1. GQ histogram ───────────────────────────────────────────────────────
    gq_arr = np.array(gq_all, dtype=float)
    ax_gq_hist.hist(gq_arr, bins=100, color="#9b59b6", alpha=0.85, edgecolor="none")
    med = np.median(gq_arr)
    ax_gq_hist.axvline(med, color=ACCENT, lw=1.8, linestyle="--",
                       label=f"median = {med:.1f}")
    ax_gq_hist.axvline(20, color="#e74c3c", lw=1.2, linestyle=":",
                       label="GQ = 20 (common filter)")
    ax_gq_hist.legend(fontsize=8, labelcolor=TEXT, facecolor=BG, edgecolor="none")
    ax_gq_hist.set_title("Genotype Quality (GQ) — all samples", fontsize=11, fontweight="bold")
    ax_gq_hist.set_xlabel("GQ", fontsize=9)
    ax_gq_hist.set_ylabel("Genotype count", fontsize=9)

    # ── 2. GQ violin per sample ───────────────────────────────────────────────
    if n_samples <= max_samples:
        rows = [(s, gq) for s, gqs in sample_gqs.items() for gq in gqs]
        df_gq = pd.DataFrame(rows, columns=["sample", "GQ"])
        order = (df_gq.groupby("sample")["GQ"]
                      .median()
                      .sort_values()
                      .index.tolist())
        sns.violinplot(data=df_gq, x="sample", y="GQ", order=order,
                       ax=ax_gq_violin, inner="quartile",
                       color="#3498db", linewidth=0.6, cut=0)
        ax_gq_violin.set_xticklabels(
            ax_gq_violin.get_xticklabels(),
            rotation=90, fontsize=max(5, 8 - n_samples // 20)
        )
        ax_gq_violin.set_title("GQ distribution per sample", fontsize=11, fontweight="bold")
        ax_gq_violin.set_xlabel("Sample", fontsize=9)
        ax_gq_violin.set_ylabel("GQ", fontsize=9)
    else:
        ax_gq_violin.text(0.5, 0.5,
                          f"Per-sample GQ violin skipped\n({n_samples} samples > --max-samples {max_samples})\n"
                          f"Re-run with --max-samples {n_samples} to enable",
                          ha="center", va="center", color="#888899",
                          fontsize=10, transform=ax_gq_violin.transAxes, linespacing=1.8)
        ax_gq_violin.set_title("GQ distribution per sample", fontsize=11, fontweight="bold")

    # ── 3. Variant count per sample — violin + strip ──────────────────────────
    df_cnt = counts_series.reset_index()
    df_cnt.columns = ["sample", "n_variants"]

    parts = ax_cnt_violin.violinplot(
        df_cnt["n_variants"], vert=False, showmedians=True,
        showextrema=True, widths=0.7
    )
    for pc in parts["bodies"]:
        pc.set_facecolor("#2ecc71")
        pc.set_alpha(0.6)
        pc.set_edgecolor("none")
    parts["cmedians"].set_color(ACCENT)
    parts["cmedians"].set_linewidth(2)
    for key in ("cmins", "cmaxes", "cbars"):
        parts[key].set_color("#888899")
        parts[key].set_linewidth(1)

    plot_samples = min(n_samples, max_samples)
    sample_subset = df_cnt.sample(n=plot_samples, random_state=42) \
                    if n_samples > plot_samples else df_cnt
    jitter = np.random.default_rng(0).uniform(-0.25, 0.25, len(sample_subset))
    ax_cnt_violin.scatter(
        sample_subset["n_variants"], 1 + jitter,
        alpha=0.5, s=18, color="#3498db", zorder=3,
        label=f"samples (n={plot_samples}{'/' + str(n_samples) if n_samples > plot_samples else ''})"
    )

    med_cnt  = df_cnt["n_variants"].median()
    mean_cnt = df_cnt["n_variants"].mean()
    ax_cnt_violin.axvline(med_cnt,  color=ACCENT,    lw=1.5, linestyle="--",
                           label=f"median = {med_cnt:,.0f}")
    ax_cnt_violin.axvline(mean_cnt, color="#e67e22", lw=1.5, linestyle=":",
                           label=f"mean = {mean_cnt:,.0f}")
    ax_cnt_violin.legend(fontsize=9, labelcolor=TEXT, facecolor=BG,
                          edgecolor="none", loc="upper right")
    ax_cnt_violin.set_title(
        f"Variant calls per sample  (n={n_samples} samples, {n_variants:,} total records)",
        fontsize=11, fontweight="bold"
    )
    ax_cnt_violin.set_xlabel("Number of variant calls (non-HOM_REF genotypes)", fontsize=9)
    ax_cnt_violin.set_yticks([])
    ax_cnt_violin.set_xlim(left=0)

    vcf_name = Path(vcf_path).name
    fig.suptitle(f"GLnexus / DeepVariant VCF Quality — {vcf_name}",
                 fontsize=13, fontweight="bold", color=TEXT, y=0.98)

    plt.savefig(out_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[✓] Plot saved → {out_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot GQ and per-sample variant counts from a GLnexus VCF."
    )
    parser.add_argument("vcf", help="Input VCF or VCF.gz")
    parser.add_argument("-o", "--output-prefix", default=None,
                        help="Output file prefix (default: derived from input)")
    parser.add_argument("--min-gq", type=float, default=0.0,
                        help="Minimum GQ to include in plots (default: 0)")
    parser.add_argument("--pass-only", action="store_true",
                        help="Only include PASS records")
    parser.add_argument("--max-samples", type=int, default=356,
                        help="Max samples for per-sample violin (default: 356)")
    args = parser.parse_args()

    prefix = args.output_prefix or Path(args.vcf).stem.replace(".vcf", "")

    gq_all, sample_counts, sample_gqs, n_variants = parse_vcf(
        args.vcf, pass_only=args.pass_only, min_gq=args.min_gq
    )

    plot_report(
        gq_all, sample_counts, sample_gqs, n_variants,
        out_path=f"{prefix}_gq_report.png",
        vcf_path=args.vcf,
        max_samples=args.max_samples,
    )


if __name__ == "__main__":
    main()
