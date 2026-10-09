#!/usr/bin/env python3
"""
compare_sv_af.py

Compare allele frequencies of matching structural variants between two VCFs,
and plot the differences.

Matching: by default, matches sites using the INFO/MatchId field (e.g.
MatchId=40.1.0) written by `truvari bench` into its tp-base.vcf.gz and
tp-comp.vcf.gz outputs - the same MatchId value appears on both sides of a
matched pair, so this is an exact join, not a distance-based guess. Run this
script directly on truvari's tp-base.vcf.gz (--vcf1) and tp-comp.vcf.gz
(--vcf2) outputs.

A --match-mode position fallback is also available for VCFs that don't have
MatchId (e.g. if you haven't run truvari, or are comparing two VCFs some
other way) - see --pos-window / --ignore-svtype below. Position matching is
approximate; MatchId matching is exact and is the recommended default.

Allele frequency, for each site, is:
  - af_info      : taken directly from INFO/AF, if present (NaN otherwise)
  - af_computed  : (alt allele count across all genotypes) / (2 x N samples),
                    counting hom-ref (0/0) and missing (./.) genotypes in the
                    denominator.
By default the plots use af_info when present, falling back to af_computed
per-site if INFO/AF is missing for that record. Use --af-source to force one
or the other explicitly.

Outputs (written to --outdir):
  - matched_af.csv           : one row per matched SV pair, with af1, af2, delta
  - af_comparison_scatter.png: AF (vcf1) vs AF (vcf2) scatter, y=x reference line
  - af_delta_histogram.png   : histogram of af2 - af1
  - af_delta_by_svtype.png   : delta AF distribution, faceted by SV type
  - summary_stats.txt        : match counts, correlation, mean/median delta

Usage:
    # recommended: run on truvari bench's own output pair
    python compare_sv_af.py --vcf1 result/tp-base.vcf.gz --vcf2 result/tp-comp.vcf.gz \\
        --outdir af_compare_out --label1 Baseline --label2 Comparison

    # fallback, no MatchId available
    python compare_sv_af.py --vcf1 a.vcf.gz --vcf2 b.vcf.gz --outdir out \\
        --match-mode position --pos-window 50 --af-source computed
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

try:
    import cyvcf2
except ImportError:
    sys.exit("ERROR: cyvcf2 is required. Install with:\n  pip install cyvcf2 --break-system-packages")

sns.set_theme(style="whitegrid")


def parse_args():
    p = argparse.ArgumentParser(description="Compare AF of matching SVs between two VCFs.")
    p.add_argument("--vcf1", required=True, help="First VCF (e.g. truvari bench's tp-base.vcf.gz).")
    p.add_argument("--vcf2", required=True, help="Second VCF (e.g. truvari bench's tp-comp.vcf.gz).")
    p.add_argument("--label1", default=None, help="Display name for VCF #1 (default: filename).")
    p.add_argument("--label2", default=None, help="Display name for VCF #2 (default: filename).")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--match-mode", choices=["matchid", "position"], default="matchid",
                    help="matchid (default): join on INFO/MatchId (e.g. from truvari bench). "
                         "position: fall back to chrom/pos(+svtype) proximity matching.")
    p.add_argument("--match-field", default="MatchId",
                    help="INFO field to join on when --match-mode matchid (default: MatchId).")
    p.add_argument("--pos-window", type=int, default=0,
                    help="[position mode only] Max distance (bp) between POS to call a match (default: 0, exact).")
    p.add_argument("--require-same-svtype", dest="require_svtype", action="store_true", default=True,
                    help="[position mode only] Require SVTYPE to match too (default: on).")
    p.add_argument("--ignore-svtype", dest="require_svtype", action="store_false",
                    help="[position mode only] Match on position alone, ignoring SVTYPE.")
    p.add_argument("--af-source", choices=["auto", "info", "computed"], default="auto",
                    help="auto (default): use INFO/AF, fall back to computed per-site if missing. "
                         "info: only use INFO/AF (unmatched sites without it are dropped). "
                         "computed: always use the counted alt-allele/2N definition.")
    return p.parse_args()


def load_sites(vcf_path, match_field=None):
    """Return a DataFrame with one row per SV site: chrom, pos, svtype, af_info, af_computed,
    and (if match_field given) that INFO field's value as 'match_id'."""
    vcf = cyvcf2.VCF(vcf_path)
    n_samples = len(vcf.samples)
    rows = []
    for v in vcf:
        svtype = v.INFO.get("SVTYPE")
        if svtype is None:
            alt = v.ALT[0] if v.ALT else ""
            svtype = alt.strip("<>") if alt.startswith("<") else "UNKNOWN"

        af_info = v.INFO.get("AF")
        if isinstance(af_info, (list, tuple)):
            af_info = af_info[0]

        af_computed = None
        if n_samples > 0:
            gt_types = v.gt_types
            n_alt = sum(1 if gt == 1 else 2 if gt == 3 else 0 for gt in gt_types)
            af_computed = n_alt / (2 * n_samples)

        row = {
            "chrom": v.CHROM,
            "pos": v.POS,
            "svtype": svtype,
            "af_info": af_info,
            "af_computed": af_computed,
        }
        if match_field:
            row["match_id"] = v.INFO.get(match_field)
        rows.append(row)
    df = pd.DataFrame(rows)
    return df


def resolve_af(df, source):
    if source == "info":
        df["af"] = df["af_info"]
    elif source == "computed":
        df["af"] = df["af_computed"]
    else:  # auto
        df["af"] = df["af_info"].where(df["af_info"].notna(), df["af_computed"])
    return df


def match_by_id(df1, df2):
    """Join df1 and df2 on their match_id column (exact match, e.g. truvari's MatchId)."""
    n1_missing = df1["match_id"].isna().sum()
    n2_missing = df2["match_id"].isna().sum()
    if n1_missing:
        print(f"Note: {n1_missing} site(s) in VCF #1 have no match_id and will be dropped.")
    if n2_missing:
        print(f"Note: {n2_missing} site(s) in VCF #2 have no match_id and will be dropped.")

    d1 = df1.dropna(subset=["match_id"])
    d2 = df2.dropna(subset=["match_id"])

    merged = d1.merge(
        d2, on="match_id", how="inner", suffixes=("1", "2"), validate="one_to_one"
    )
    matched = pd.DataFrame({
        "match_id": merged["match_id"],
        "chrom": merged["chrom1"],
        "pos1": merged["pos1"],
        "pos2": merged["pos2"],
        "svtype": merged["svtype1"],
        "af1": merged["af1"],
        "af2": merged["af2"],
    })
    return matched


def match_by_position(df1, df2, pos_window, require_svtype):
    """Match sites in df1 to df2 by chrom (+ svtype) and nearest position within pos_window.
    Returns a DataFrame with one row per matched pair."""
    matches = []
    group_cols = ["chrom", "svtype"] if require_svtype else ["chrom"]

    for key, g1 in df1.groupby(group_cols):
        key_tuple = key if isinstance(key, tuple) else (key,)
        mask = pd.Series(True, index=df2.index)
        for col, val in zip(group_cols, key_tuple):
            mask &= (df2[col] == val)
        g2 = df2[mask]
        if g2.empty:
            continue

        g2_positions = g2["pos"].to_numpy()
        used_g2_idx = set()

        for _, row1 in g1.sort_values("pos").iterrows():
            if len(used_g2_idx) == len(g2):
                break
            diffs = np.abs(g2_positions - row1["pos"])
            diffs_masked = diffs.copy().astype(float)
            for used_pos_idx in used_g2_idx:
                diffs_masked[used_pos_idx] = np.inf
            best_idx = int(np.argmin(diffs_masked))
            best_dist = diffs_masked[best_idx]
            if best_dist <= pos_window:
                row2 = g2.iloc[best_idx]
                used_g2_idx.add(best_idx)
                matches.append({
                    "chrom": row1["chrom"],
                    "pos1": row1["pos"],
                    "pos2": row2["pos"],
                    "svtype": row1["svtype"] if require_svtype else f"{row1['svtype']}/{row2['svtype']}",
                    "af1": row1["af"],
                    "af2": row2["af"],
                })

    return pd.DataFrame(matches)


def make_plots(matched_df, label1, label2, outdir):
    os.makedirs(outdir, exist_ok=True)
    df = matched_df.dropna(subset=["af1", "af2"]).copy()
    df["delta"] = df["af2"] - df["af1"]

    if df.empty:
        print("No matched sites have usable AF values on both sides - skipping plots.")
        return df

    # ---- 1. scatter: af1 vs af2, y=x reference ----
    plt.figure(figsize=(6.5, 6.5))
    ax = sns.scatterplot(data=df, x="af1", y="af2", hue="svtype", alpha=0.6, s=25)
    ax.plot([0, 1], [0, 1], color="black", linestyle="--", linewidth=1, label="y = x")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(f"AF - {label1}")
    ax.set_ylabel(f"AF - {label2}")
    ax.set_title("Matched SV allele frequency comparison")
    ax.legend(title="SVTYPE", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "af_comparison_scatter.png"), dpi=150)
    plt.close()

    # ---- 2. delta AF histogram (all types combined) ----
    plt.figure(figsize=(7, 5))
    ax = sns.histplot(data=df, x="delta", bins=40, color="teal")
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    mean_delta = df["delta"].mean()
    ax.axvline(mean_delta, color="darkorange", linestyle="-", linewidth=1.5,
               label=f"mean = {mean_delta:.3f}")
    ax.set_xlabel(f"AF({label2}) - AF({label1})")
    ax.set_ylabel("Number of matched SVs")
    ax.set_title("Distribution of AF differences")
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "af_delta_histogram.png"), dpi=150)
    plt.close()

    # ---- 3. delta AF by SVTYPE ----
    if df["svtype"].nunique() > 1:
        plt.figure(figsize=(max(7, df["svtype"].nunique() * 1.3), 5))
        order = sorted(df["svtype"].unique())
        ax = sns.violinplot(data=df, x="svtype", y="delta", order=order, hue="svtype",
                             palette="Set2", legend=False, cut=0)
        sns.stripplot(data=df, x="svtype", y="delta", order=order, color="black", alpha=0.3, size=2, ax=ax)
        ax.axhline(0, color="black", linestyle="--", linewidth=1)
        ax.set_xlabel("SVTYPE")
        ax.set_ylabel(f"AF({label2}) - AF({label1})")
        ax.set_title("AF difference by SV type")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "af_delta_by_svtype.png"), dpi=150)
        plt.close()

    df.to_csv(os.path.join(outdir, "matched_af.csv"), index=False)
    return df


def write_summary(df1_n, df2_n, matched_df, df, label1, label2, outdir, match_mode):
    lines = []
    lines.append("SV allele frequency comparison summary")
    lines.append("=" * 45)
    lines.append(f"{label1}: {df1_n} total SV sites")
    lines.append(f"{label2}: {df2_n} total SV sites")
    if match_mode == "matchid":
        lines.append(f"Matched pairs (joined on MatchId): {len(matched_df)}")
    else:
        lines.append(f"Matched pairs (position +/- window, same chrom/svtype): {len(matched_df)}")
    lines.append(f"Matched pairs with usable AF on both sides: {len(df)}")
    lines.append("")

    if not df.empty:
        corr = df["af1"].corr(df["af2"])
        spearman = df["af1"].corr(df["af2"], method="spearman")
        lines.append(f"Pearson correlation (af1 vs af2):  {corr:.4f}")
        lines.append(f"Spearman correlation (af1 vs af2): {spearman:.4f}")
        lines.append("")
        lines.append("AF delta (af2 - af1) stats:")
        lines.append(f"  mean:   {df['delta'].mean():.4f}")
        lines.append(f"  median: {df['delta'].median():.4f}")
        lines.append(f"  std:    {df['delta'].std():.4f}")
        lines.append(f"  min:    {df['delta'].min():.4f}")
        lines.append(f"  max:    {df['delta'].max():.4f}")
        lines.append("")
        lines.append("Mean |delta| by SVTYPE:")
        for svtype, sub in df.groupby("svtype"):
            lines.append(f"  {svtype:10s} n={len(sub):5d}  mean|delta|={sub['delta'].abs().mean():.4f}")

    text = "\n".join(lines)
    with open(os.path.join(outdir, "summary_stats.txt"), "w") as f:
        f.write(text + "\n")
    print(text)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    label1 = args.label1 or os.path.basename(args.vcf1)
    label2 = args.label2 or os.path.basename(args.vcf2)

    match_field = args.match_field if args.match_mode == "matchid" else None

    print(f"Loading {args.vcf1} ...")
    df1 = load_sites(args.vcf1, match_field=match_field)
    print(f"  {len(df1)} sites")

    print(f"Loading {args.vcf2} ...")
    df2 = load_sites(args.vcf2, match_field=match_field)
    print(f"  {len(df2)} sites")

    df1 = resolve_af(df1, args.af_source)
    df2 = resolve_af(df2, args.af_source)

    if args.match_mode == "matchid":
        print(f"Matching sites by INFO/{args.match_field} ...")
        if df1["match_id"].isna().all() or df2["match_id"].isna().all():
            sys.exit(
                f"ERROR: no records in one or both VCFs have INFO/{args.match_field} set.\n"
                f"This field is written by `truvari bench` into its tp-base.vcf.gz / tp-comp.vcf.gz "
                f"outputs - point --vcf1/--vcf2 at those files, or use --match-mode position instead."
            )
        matched_df = match_by_id(df1, df2)
    else:
        print("Matching sites by position" + (" and SVTYPE" if args.require_svtype else "") +
              f" (window = {args.pos_window} bp) ...")
        matched_df = match_by_position(df1, df2, args.pos_window, args.require_svtype)

    print(f"  {len(matched_df)} matched pairs")

    if matched_df.empty:
        sys.exit("No matches found. In matchid mode, check that both VCFs actually share MatchId "
                  "values (e.g. they're the tp-base/tp-comp pair from the same truvari bench run). "
                  "In position mode, check chrom naming (e.g. 'chr1' vs '1') and --pos-window.")

    df = make_plots(matched_df, label1, label2, args.outdir)
    write_summary(len(df1), len(df2), matched_df, df, label1, label2, args.outdir, args.match_mode)

    print(f"\nDone. Outputs written to: {args.outdir}")


if __name__ == "__main__":
    main()