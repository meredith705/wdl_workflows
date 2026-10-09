#!/usr/bin/env python3
"""
plot_sv_stats.py

Summarize and plot structural variant (SV) statistics from a multisample VCF.

Expects standard SV caller INFO fields (Manta/DELLY/smoove/lumpy-style):
  - SVTYPE  (DEL, DUP, INS, INV, BND, ...)
  - SVLEN   (signed length; optional, falls back to END - POS)
  - END     (end position; optional)

Outputs (written to --outdir):
  - sv_type_counts.csv           : count of each SVTYPE (unique sites, not per-genotype)
  - sv_per_sample.csv            : number of non-ref SV genotypes per sample
  - summary_stats.txt            : plain-text overview (totals, medians, etc.)
  - sv_type_frequency.png        : bar chart of SV type counts
  - sv_length_distribution.png   : histogram of SV lengths, faceted by type
  - sv_per_sample.png            : bar chart of SV count per sample
  - sv_type_by_sample_heatmap.png: heatmap of SV type counts per sample

Usage:
    python plot_sv_stats.py --vcf cohort.sv.vcf.gz --outdir sv_stats_out
    python plot_sv_stats.py --vcf cohort.sv.vcf.gz --outdir sv_stats_out --min-len 50 --max-len 1000000
"""

import argparse
import os
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

try:
    import cyvcf2
except ImportError:
    sys.exit(
        "ERROR: cyvcf2 is required. Install with:\n"
        "  pip install cyvcf2 --break-system-packages"
    )

sns.set_theme(style="whitegrid")


def parse_args():
    p = argparse.ArgumentParser(description="Plot SV frequency and sample-level stats from a multisample VCF.")
    p.add_argument("--vcf", required=True, help="Path to multisample SV VCF (bgzipped or plain).")
    p.add_argument("--outdir", required=True, help="Directory to write plots and tables to.")
    p.add_argument("--min-len", type=int, default=None,
                    help="Exclude SVs with |SVLEN| below this (bp). BND (no length) always kept.")
    p.add_argument("--max-len", type=int, default=None,
                    help="Exclude SVs with |SVLEN| above this (bp). BND (no length) always kept.")
    p.add_argument("--pass-only", action="store_true",
                    help="Only include variants with FILTER == PASS (or missing/'.').")
    return p.parse_args()


def load_variants(vcf_path, min_len=None, max_len=None, pass_only=False):
    """Stream the VCF once, returning:
       - a per-site dataframe (site_id, chrom, pos, svtype, svlen)
       - a per-sample non-ref SV-genotype count dict
       - a per-sample x svtype count dataframe
    """
    vcf = cyvcf2.VCF(vcf_path)
    samples = vcf.samples

    site_rows = []
    sample_counts = {s: 0 for s in samples}
    sample_type_counts = {s: {} for s in samples}

    n_total = 0
    n_kept = 0

    for v in vcf:
        n_total += 1

        if pass_only and v.FILTER not in (None, "PASS", "."):
            continue

        svtype = v.INFO.get("SVTYPE")
        if svtype is None:
            # fall back to inferring from ALT allele, e.g. "<DEL>"
            alt = v.ALT[0] if v.ALT else ""
            svtype = alt.strip("<>") if alt.startswith("<") else "UNKNOWN"

        svlen = v.INFO.get("SVLEN")
        if svlen is None:
            end = v.INFO.get("END")
            svlen = (end - v.POS) if end is not None else None
        if isinstance(svlen, (list, tuple)):
            svlen = svlen[0]

        abs_len = abs(svlen) if svlen is not None else None

        if abs_len is not None:
            if min_len is not None and abs_len < min_len:
                continue
            if max_len is not None and abs_len > max_len:
                continue

        n_kept += 1
        site_rows.append({
            "site_id": v.ID if v.ID else f"{v.CHROM}:{v.POS}",
            "chrom": v.CHROM,
            "pos": v.POS,
            "svtype": svtype,
            "svlen": abs_len,
        })

        # per-sample genotype counting: count non-ref, non-missing calls
        gt_types = v.gt_types  # 0=HOM_REF, 1=HET, 2=UNKNOWN/missing(varies by cyvcf2 version), 3=HOM_ALT
        for s, gt in zip(samples, gt_types):
            if gt in (1, 3):  # het or hom-alt -> sample carries this SV
                sample_counts[s] += 1
                sample_type_counts[s][svtype] = sample_type_counts[s].get(svtype, 0) + 1

    print(f"Parsed {n_total} records, kept {n_kept} after filters.")

    site_df = pd.DataFrame(site_rows)
    sample_df = pd.DataFrame({
        "sample": list(sample_counts.keys()),
        "n_svs": list(sample_counts.values()),
    }).sort_values("n_svs", ascending=False).reset_index(drop=True)

    type_by_sample_df = pd.DataFrame(sample_type_counts).T.fillna(0).astype(int)
    type_by_sample_df.index.name = "sample"

    return site_df, sample_df, type_by_sample_df


def make_plots(site_df, sample_df, type_by_sample_df, outdir):
    os.makedirs(outdir, exist_ok=True)

    # ---- 1. SV type frequency (unique sites) ----
    type_counts = site_df["svtype"].value_counts().sort_values(ascending=False)
    type_counts.to_csv(os.path.join(outdir, "sv_type_counts.csv"), header=["count"])

    plt.figure(figsize=(7, 5))
    ax = sns.barplot(x=type_counts.index, y=type_counts.values, hue=type_counts.index,
                      palette="viridis", legend=False)
    ax.set_xlabel("SV type")
    ax.set_ylabel("Number of sites")
    ax.set_title("SV type frequency (unique sites)")
    for i, v in enumerate(type_counts.values):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "sv_type_frequency.png"), dpi=150)
    plt.close()

    # ---- 2. SV length distribution, per type (excludes BND / length-less types) ----
    len_df = site_df.dropna(subset=["svlen"])
    if not len_df.empty:
        g = sns.displot(
            data=len_df, x="svlen", col="svtype", col_wrap=3,
            bins=40, log_scale=(True, False), facet_kws={"sharex": False, "sharey": False},
            color="steelblue",
        )
        g.set_axis_labels("SV length (bp, log scale)", "Count")
        g.fig.suptitle("SV length distribution by type", y=1.02)
        g.savefig(os.path.join(outdir, "sv_length_distribution.png"), dpi=150, bbox_inches="tight")
        plt.close(g.fig)
    else:
        print("No SVs with a usable length (all BND or missing SVLEN/END) - skipping length plot.")

    # ---- 3. Number of SVs per sample ----
    sample_df.to_csv(os.path.join(outdir, "sv_per_sample.csv"), index=False)

    plt.figure(figsize=(max(8, len(sample_df) * 0.35), 5))
    ax = sns.barplot(data=sample_df, x="sample", y="n_svs", color="darkorange")
    ax.set_xlabel("Sample")
    ax.set_ylabel("Number of SVs (het + hom-alt genotypes)")
    ax.set_title("SVs per sample")
    plt.xticks(rotation=90)
    mean_val = sample_df["n_svs"].mean()
    ax.axhline(mean_val, color="black", linestyle="--", linewidth=1, label=f"mean = {mean_val:.1f}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "sv_per_sample.png"), dpi=150)
    plt.close()

    # ---- 4. Heatmap: SV type counts per sample ----
    if not type_by_sample_df.empty:
        plt.figure(figsize=(max(8, type_by_sample_df.shape[0] * 0.4), max(4, type_by_sample_df.shape[1] * 0.5)))
        sns.heatmap(type_by_sample_df.T, cmap="mako", annot=False, cbar_kws={"label": "SV count"})
        plt.xlabel("Sample")
        plt.ylabel("SV type")
        plt.title("SV type counts per sample")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "sv_type_by_sample_heatmap.png"), dpi=150)
        plt.close()

    return type_counts


def write_summary(site_df, sample_df, type_counts, outdir):
    lines = []
    lines.append("SV summary statistics")
    lines.append("=" * 40)
    lines.append(f"Total unique SV sites: {len(site_df)}")
    lines.append("")
    lines.append("SV type breakdown:")
    for t, c in type_counts.items():
        pct = 100 * c / len(site_df) if len(site_df) else 0
        lines.append(f"  {t:10s} {c:6d}  ({pct:.1f}%)")
    lines.append("")

    len_df = site_df.dropna(subset=["svlen"])
    if not len_df.empty:
        lines.append("SV length stats (bp, absolute value, excludes BND):")
        lines.append(f"  min:    {len_df['svlen'].min():.0f}")
        lines.append(f"  median: {len_df['svlen'].median():.0f}")
        lines.append(f"  mean:   {len_df['svlen'].mean():.0f}")
        lines.append(f"  max:    {len_df['svlen'].max():.0f}")
        lines.append("")

    lines.append("Per-sample SV counts:")
    lines.append(f"  n_samples: {len(sample_df)}")
    lines.append(f"  min:    {sample_df['n_svs'].min()}")
    lines.append(f"  median: {sample_df['n_svs'].median():.1f}")
    lines.append(f"  mean:   {sample_df['n_svs'].mean():.1f}")
    lines.append(f"  max:    {sample_df['n_svs'].max()}")
    lines.append(f"  std:    {sample_df['n_svs'].std():.1f}")
    lines.append("")
    lines.append("Top 5 samples by SV count:")
    for _, row in sample_df.head(5).iterrows():
        lines.append(f"  {row['sample']:20s} {row['n_svs']}")
    lines.append("")
    lines.append("Bottom 5 samples by SV count:")
    for _, row in sample_df.tail(5).iterrows():
        lines.append(f"  {row['sample']:20s} {row['n_svs']}")

    text = "\n".join(lines)
    with open(os.path.join(outdir, "summary_stats.txt"), "w") as f:
        f.write(text + "\n")
    print(text)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    site_df, sample_df, type_by_sample_df = load_variants(
        args.vcf, min_len=args.min_len, max_len=args.max_len, pass_only=args.pass_only
    )

    if site_df.empty:
        sys.exit("No SV sites survived parsing/filtering - nothing to plot.")

    type_counts = make_plots(site_df, sample_df, type_by_sample_df, args.outdir)
    write_summary(site_df, sample_df, type_counts, args.outdir)

    print(f"\nDone. Outputs written to: {args.outdir}")


if __name__ == "__main__":
    main()