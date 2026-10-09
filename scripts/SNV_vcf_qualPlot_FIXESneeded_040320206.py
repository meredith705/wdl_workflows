#!/usr/bin/env python3
"""
plot_vcf_quality.py — Plot quality metrics from a SNV VCF file.

Usage:
    python plot_vcf_quality.py input.vcf[.gz] [-o output_prefix] [--min-qual 0] [--pass-only]

Outputs:
    <prefix>_quality_report.png  — multi-panel figure
    <prefix>_quality_stats.tsv   — summary statistics table

Requirements:
    pip install cyvcf2 matplotlib seaborn pandas numpy
"""

import argparse
import sys
from pathlib import Path

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


# ── helpers ──────────────────────────────────────────────────────────────────

def safe_info(variant, key, default=np.nan):
    try:
        val = variant.INFO.get(key)
        return float(val) if val is not None else default
    except (TypeError, ValueError):
        return default


def parse_vcf(vcf_path: str, pass_only: bool, min_qual: float) -> pd.DataFrame:
    """Extract per-variant quality metrics into a DataFrame."""
    records = []
    vcf = VCF(vcf_path)

    for v in vcf:
        # SNVs only: ref and alt must each be a single base
        if len(v.REF) != 1:
            continue
        alts = v.ALT
        if not alts or any(len(a) != 1 or a == "*" for a in alts):
            continue

        # FILTER
        filt = v.FILTER  # None → PASS in cyvcf2
        filter_str = filt if filt else "PASS"
        if pass_only and filter_str != "PASS":
            continue

        qual = v.QUAL if v.QUAL is not None else np.nan
        if not np.isnan(qual) and qual < min_qual:
            continue

        records.append({
            "CHROM":  v.CHROM,
            "POS":    v.POS,
            "REF":    v.REF,
            "ALT":    ",".join(alts),
            "QUAL":   qual,
            "FILTER": filter_str,
            # Common INFO fields (present in GATK, DeepVariant, etc.)
            "DP":     safe_info(v, "DP"),
            "QD":     safe_info(v, "QD"),
            "MQ":     safe_info(v, "MQ"),
            "FS":     safe_info(v, "FS"),
            "SOR":    safe_info(v, "SOR"),
            "MQRankSum": safe_info(v, "MQRankSum"),
            "ReadPosRankSum": safe_info(v, "ReadPosRankSum"),
            "AF":     safe_info(v, "AF"),
        })

    vcf.close()

    if not records:
        sys.exit("No SNV records found after filtering. Check --pass-only / --min-qual.")

    df = pd.DataFrame(records)
    # Derive transition/transversion label
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    df["VARIANT_TYPE"] = df.apply(
        lambda r: "Ti" if (r.REF, r.ALT) in transitions else "Tv", axis=1
    )
    return df


# ── plotting ─────────────────────────────────────────────────────────────────

PALETTE = {"PASS": "#2ecc71", "FAIL": "#e74c3c", "Ti": "#3498db", "Tv": "#e67e22"}

def plot_report(df: pd.DataFrame, out_path: str, vcf_path: str):
    total = len(df)
    pass_n = (df["FILTER"] == "PASS").sum()
    ti = (df["VARIANT_TYPE"] == "Ti").sum()
    tv = (df["VARIANT_TYPE"] == "Tv").sum()
    ti_tv = ti / tv if tv else np.nan

    fig = plt.figure(figsize=(18, 14))
    fig.patch.set_facecolor("#1a1a2e")
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

    ax_qual   = fig.add_subplot(gs[0, 0])
    ax_dp     = fig.add_subplot(gs[0, 1])
    ax_qd     = fig.add_subplot(gs[0, 2])
    ax_mq     = fig.add_subplot(gs[1, 0])
    ax_fs     = fig.add_subplot(gs[1, 1])
    ax_titv   = fig.add_subplot(gs[1, 2])
    ax_filter = fig.add_subplot(gs[2, 0])
    ax_scatter= fig.add_subplot(gs[2, 1])
    ax_stats  = fig.add_subplot(gs[2, 2])

    axes = [ax_qual, ax_dp, ax_qd, ax_mq, ax_fs,
            ax_titv, ax_filter, ax_scatter, ax_stats]
    for ax in axes:
        ax.set_facecolor("#16213e")
        for spine in ax.spines.values():
            spine.set_color("#444466")
        ax.tick_params(colors="#ccccdd", labelsize=8)
        ax.xaxis.label.set_color("#ccccdd")
        ax.yaxis.label.set_color("#ccccdd")
        ax.title.set_color("#eeeeff")

    def _hist(ax, series, title, xlabel, color, log_scale=False):
        data = series.dropna()
        if data.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    color="#888899", transform=ax.transAxes)
            ax.set_title(title, fontsize=10, fontweight="bold")
            return
        ax.hist(data, bins=60, color=color, alpha=0.85, edgecolor="none")
        med = data.median()
        ax.axvline(med, color="#ffdd57", lw=1.5, linestyle="--",
                   label=f"median={med:.1f}")
        ax.legend(fontsize=7, labelcolor="#ffdd57",
                  facecolor="#1a1a2e", edgecolor="none")
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.set_xlabel(xlabel, fontsize=8)
        ax.set_ylabel("Count", fontsize=8)
        if log_scale:
            ax.set_yscale("log")

    # 1. QUAL
    _hist(ax_qual, df["QUAL"], "Variant QUAL Score", "QUAL", "#9b59b6")

    # 2. DP
    dp = df["DP"].dropna()
    if not dp.empty:
        cap = np.percentile(dp, 99)
        _hist(ax_dp, dp.clip(upper=cap), f"Depth (DP) [capped p99={cap:.0f}]", "DP", "#1abc9c")
    else:
        ax_dp.text(0.5, 0.5, "DP not in INFO", ha="center", va="center",
                   color="#888899", transform=ax_dp.transAxes)
        ax_dp.set_title("Depth (DP)", fontsize=10, fontweight="bold")

    # 3. QD
    _hist(ax_qd, df["QD"], "Quality by Depth (QD)", "QD", "#3498db")

    # 4. MQ
    _hist(ax_mq, df["MQ"], "Mapping Quality (MQ)", "MQ", "#e67e22")

    # 5. FS (log scale)
    _hist(ax_fs, df["FS"], "Fisher Strand Bias (FS)", "FS", "#e74c3c", log_scale=True)

    # 6. Ti/Tv bar
    titv_data = df["VARIANT_TYPE"].value_counts()
    colors_titv = [PALETTE.get(k, "#aaaaaa") for k in titv_data.index]
    bars = ax_titv.bar(titv_data.index, titv_data.values,
                       color=colors_titv, edgecolor="none", width=0.5)
    for bar, val in zip(bars, titv_data.values):
        ax_titv.text(bar.get_x() + bar.get_width() / 2,
                     bar.get_height() + total * 0.01,
                     f"{val:,}", ha="center", va="bottom",
                     color="#eeeeff", fontsize=9)
    ax_titv.set_title(f"Ti/Tv Ratio: {ti_tv:.3f}", fontsize=10, fontweight="bold")
    ax_titv.set_ylabel("Count", fontsize=8)

    # 7. FILTER breakdown
    filter_counts = df["FILTER"].value_counts()
    fc = [PALETTE.get(k, "#aaaaaa") for k in filter_counts.index]
    ax_filter.barh(filter_counts.index, filter_counts.values,
                   color=fc, edgecolor="none")
    ax_filter.set_title("FILTER Status", fontsize=10, fontweight="bold")
    ax_filter.set_xlabel("Count", fontsize=8)

    # 8. QUAL vs DP scatter
    scatter_df = df[["QUAL", "DP", "VARIANT_TYPE"]].dropna()
    if not scatter_df.empty:
        dp_cap = np.percentile(scatter_df["DP"], 99)
        scatter_df = scatter_df[scatter_df["DP"] <= dp_cap]
        for vt, grp in scatter_df.groupby("VARIANT_TYPE"):
            ax_scatter.scatter(grp["DP"], grp["QUAL"],
                               alpha=0.15, s=4,
                               color=PALETTE.get(vt, "#aaaaaa"),
                               label=vt, rasterized=True)
        ax_scatter.set_xlabel("DP", fontsize=8)
        ax_scatter.set_ylabel("QUAL", fontsize=8)
        ax_scatter.set_title("QUAL vs DP", fontsize=10, fontweight="bold")
        ax_scatter.legend(fontsize=7, facecolor="#1a1a2e",
                          edgecolor="none", labelcolor="#eeeeff")
    else:
        ax_scatter.text(0.5, 0.5, "Insufficient data",
                        ha="center", va="center", color="#888899",
                        transform=ax_scatter.transAxes)
        ax_scatter.set_title("QUAL vs DP", fontsize=10, fontweight="bold")

    # 9. Summary stats table
    ax_stats.axis("off")
    metrics = {
        "Total SNVs":    f"{total:,}",
        "PASS":          f"{pass_n:,} ({100*pass_n/total:.1f}%)",
        "Ti":            f"{ti:,}",
        "Tv":            f"{tv:,}",
        "Ti/Tv":         f"{ti_tv:.3f}" if not np.isnan(ti_tv) else "N/A",
        "Median QUAL":   f"{df['QUAL'].median():.1f}" if df['QUAL'].notna().any() else "N/A",
        "Median DP":     f"{df['DP'].median():.1f}"   if df['DP'].notna().any()   else "N/A",
        "Median QD":     f"{df['QD'].median():.1f}"   if df['QD'].notna().any()   else "N/A",
        "Median MQ":     f"{df['MQ'].median():.1f}"   if df['MQ'].notna().any()   else "N/A",
    }
    rows = list(metrics.items())
    tbl = ax_stats.table(
        cellText=rows,
        colLabels=["Metric", "Value"],
        cellLoc="left",
        loc="center",
        bbox=[0, 0, 1, 1],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_facecolor("#16213e" if r % 2 == 0 else "#1e2a4a")
        cell.set_edgecolor("#333355")
        cell.set_text_props(color="#eeeeff")
        if r == 0:
            cell.set_facecolor("#2c2c5e")
            cell.set_text_props(color="#ffdd57", fontweight="bold")
    ax_stats.set_title("Summary Statistics", fontsize=10,
                        fontweight="bold", color="#eeeeff")

    # Super-title
    vcf_name = Path(vcf_path).name
    fig.suptitle(f"SNV Quality Report — {vcf_name}",
                 fontsize=14, fontweight="bold",
                 color="#eeeeff", y=0.98)

    plt.savefig(out_path, dpi=150, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[✓] Plot saved → {out_path}")


def save_stats(df: pd.DataFrame, out_path: str):
    cols = ["QUAL", "DP", "QD", "MQ", "FS", "SOR", "MQRankSum", "ReadPosRankSum"]
    present = [c for c in cols if df[c].notna().any()]
    stats = df[present].describe().T
    stats.to_csv(out_path, sep="\t", float_format="%.4f")
    print(f"[✓] Stats saved → {out_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot quality metrics from a SNV VCF."
    )
    parser.add_argument("vcf", help="Input VCF or VCF.gz")
    parser.add_argument("-o", "--output-prefix", default=None,
                        help="Output file prefix (default: derived from input)")
    parser.add_argument("--min-qual", type=float, default=0.0,
                        help="Minimum QUAL score to include (default: 0)")
    parser.add_argument("--pass-only", action="store_true",
                        help="Only include PASS variants")
    args = parser.parse_args()

    prefix = args.output_prefix or Path(args.vcf).stem.replace(".vcf", "")

    print(f"[…] Parsing {args.vcf}")
    df = parse_vcf(args.vcf, pass_only=args.pass_only, min_qual=args.min_qual)
    print(f"[…] {len(df):,} SNVs loaded")

    plot_report(df, f"{prefix}_quality_report.png", args.vcf)
    save_stats(df,  f"{prefix}_quality_stats.tsv")


if __name__ == "__main__":
    main()