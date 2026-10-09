#!/usr/bin/env python3
"""
plot_sv_stats.py

Summarize and plot structural variant (SV) statistics from a multisample VCF.

Expects standard SV caller INFO fields (Manta/DELLY/smoove/lumpy-style):
  - SVTYPE  (DEL, DUP, INS, INV, BND, ...)
  - SVLEN   (signed length; optional, falls back to END - POS)
  - END     (end position; optional)
  - AF      (allele frequency; optional - see below for how this is used)

Two allele frequencies are computed per site, plotted separately:
  - af_vcf      : taken directly from INFO/AF, as the caller reported it (NaN if absent)
  - af_computed : (count of alt alleles across all genotypes) / (2 x number of samples in the VCF).
                  Every sample counts toward the denominator, including hom-ref (0/0) and
                  missing (./.) genotypes - missing calls simply contribute 0 alt alleles to the
                  numerator but each still counts for 2 in the denominator (diploid).

Outputs (written to --outdir):
  - sv_type_counts.csv              : count of each SVTYPE (unique sites, not per-genotype)
  - sv_per_sample.csv               : number of non-ref SV genotypes per sample
  - sv_allele_frequencies.csv       : per-site af_vcf and af_computed, alongside svtype/svlen
  - summary_stats.txt               : plain-text overview (totals, medians, etc.)
  - sv_type_frequency.png           : bar chart of SV type counts
  - sv_length_distribution.png      : histogram of SV lengths, faceted by type
  - sv_per_sample.png               : bar chart of SV count per sample
  - sv_type_by_sample_heatmap.png   : heatmap of SV type counts per sample
  - sv_af_distribution_vcf.png      : overall AF histogram, using INFO/AF
  - sv_af_by_type_vcf.png           : AF histogram faceted by SV type, using INFO/AF
  - sv_af_distribution_computed.png : overall AF histogram, using the counted-up definition above
  - sv_af_by_type_computed.png      : AF histogram faceted by SV type, using the counted-up definition

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

        # per-sample genotype counting: count non-ref, non-missing calls
        gt_types = v.gt_types  # 0=HOM_REF, 1=HET, 2=UNKNOWN/missing(varies by cyvcf2 version), 3=HOM_ALT

        # --- af_vcf: exactly what INFO/AF says, no fallback ---
        af_vcf = v.INFO.get("AF")
        if isinstance(af_vcf, (list, tuple)):
            af_vcf = af_vcf[0]  # first ALT allele's AF, in case of a leftover multiallelic record

        # --- af_computed: alt allele count / number of samples (0/0 and ./. both count in the denominator) ---
        n_alt_alleles = 0
        for gt in gt_types:
            if gt == 1:      # het
                n_alt_alleles += 1
            elif gt == 3:    # hom-alt
                n_alt_alleles += 2
            # gt == 0 (hom-ref) and gt == 2 (missing) both contribute 0 to the numerator,
            # but every sample is still counted in the denominator below
        af_computed = (n_alt_alleles / (2 * len(samples))) if len(samples) > 0 else None

        n_kept += 1
        site_rows.append({
            "site_id": v.ID if v.ID else f"{v.CHROM}:{v.POS}",
            "chrom": v.CHROM,
            "pos": v.POS,
            "svtype": svtype,
            "svlen": abs_len,
            "af_vcf": af_vcf,
            "af_computed": af_computed,
        })

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

    # ---- 2b. Allele frequency distributions (two versions: INFO/AF, and counted-up) ----
    af_cols = [c for c in ("af_vcf", "af_computed") if c in site_df.columns]
    combined_af_df = site_df[["svtype"] + af_cols].copy()
    combined_af_df.to_csv(os.path.join(outdir, "sv_allele_frequencies.csv"), index=False)

    def plot_af(col, label, tag, color):
        af_df = site_df.dropna(subset=[col])
        if af_df.empty:
            print(f"No SVs with a usable {col} - skipping {tag} AF plots.")
            return

        x_max = max(1.0, af_df[col].max())

        # overall AF histogram, all SV types combined, log-scaled y-axis (counts)
        plt.figure(figsize=(7, 5))
        ax = sns.histplot(data=af_df, x=col, bins=40, color=color)
        ax.set_xlabel(label)
        ax.set_ylabel("Number of SV sites (log scale)")
        ax.set_title(f"SV allele frequency distribution (all types) - {tag}")
        ax.set_xlim(0, x_max)
        ax.set_yscale("log")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"sv_af_distribution_{tag}.png"), dpi=150)
        plt.close()

        # AF histogram, faceted by SV type, log-scaled y-axis (counts)
        g = sns.displot(
            data=af_df, x=col, col="svtype", col_wrap=3,
            bins=30, facet_kws={"sharex": True, "sharey": False},
            color=color,
        )
        g.set_axis_labels(label, "Count (log scale)")
        g.set(xlim=(0, x_max))
        for ax in g.axes.flat:
            ax.set_yscale("log")
        g.fig.suptitle(f"SV allele frequency distribution by type - {tag}", y=1.02)
        g.savefig(os.path.join(outdir, f"sv_af_by_type_{tag}.png"), dpi=150, bbox_inches="tight")
        plt.close(g.fig)

    if "af_vcf" in site_df.columns:
        plot_af("af_vcf", "Allele frequency (INFO/AF)", "vcf", "mediumpurple")
    if "af_computed" in site_df.columns:
        plot_af("af_computed", "Allele frequency (counted: alt alleles / 2N samples)", "computed", "darkslateblue")

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

    for col, label in [("af_vcf", "INFO/AF"), ("af_computed", "counted (alt alleles / 2N samples)")]:
        if col not in site_df.columns:
            continue
        af_df = site_df.dropna(subset=[col])
        if not af_df.empty:
            lines.append(f"Allele frequency stats [{label}] ({len(af_df)}/{len(site_df)} sites usable):")
            lines.append(f"  min:    {af_df[col].min():.4f}")
            lines.append(f"  median: {af_df[col].median():.4f}")
            lines.append(f"  mean:   {af_df[col].mean():.4f}")
            lines.append(f"  max:    {af_df[col].max():.4f}")
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