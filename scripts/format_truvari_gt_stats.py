#!/usr/bin/env python3
"""
format_truvari_gt_stats.py

Format the single-sample genotype-comparison stats from a `truvari bench`
run into a readable report, plus a heatmap of the base-vs-comparison
genotype confusion matrix (gt_matrix).

Reads summary.json as written by `truvari bench -b base.vcf.gz -c comp.vcf.gz
-o result/` (this script accepts either --summary result/summary.json
directly, or --bench-dir result/ to find it automatically).

Outputs (written to --outdir):
  - gt_stats_report.txt  : plain-text formatted report (also printed to stdout)
  - gt_matrix.csv         : the base x comp genotype confusion matrix as a table
  - gt_matrix_heatmap.png : heatmap of the confusion matrix (counts + row-normalized %)

Usage:
    python format_truvari_gt_stats.py --bench-dir result/ --outdir gt_report
    python format_truvari_gt_stats.py --summary result/summary.json --outdir gt_report
"""

import argparse
import json
import os
import sys

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="white")


def parse_args():
    p = argparse.ArgumentParser(description="Format truvari bench single-sample GT comparison stats.")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--summary", help="Path directly to a truvari bench summary.json.")
    g.add_argument("--bench-dir", help="Path to a truvari bench output directory (looks for summary.json inside).")
    p.add_argument("--outdir", required=True, help="Directory to write the report/plot to.")
    p.add_argument("--label", default=None, help="Optional label for the report title (e.g. sample name).")
    return p.parse_args()


def load_summary(args):
    if args.summary:
        path = args.summary
    else:
        path = os.path.join(args.bench_dir, "summary.json")
    if not os.path.isfile(path):
        sys.exit(f"ERROR: summary.json not found at {path}")
    with open(path) as f:
        return json.load(f), path


def fmt_pct(x):
    return f"{x * 100:.2f}%" if x is not None else "n/a"


def fmt_num(x):
    return f"{x:,}" if x is not None else "n/a"


def build_gt_matrix_df(gt_matrix):
    """gt_matrix in summary.json is {base_gt: {comp_gt: count}}. Build a full, symmetric,
    zero-filled DataFrame (rows=base GT, cols=comp GT) so nothing is silently missing."""
    if not gt_matrix:
        return None
    all_gts = sorted(set(gt_matrix.keys()) | {c for row in gt_matrix.values() for c in row.keys()})
    df = pd.DataFrame(0, index=all_gts, columns=all_gts, dtype=int)
    for base_gt, row in gt_matrix.items():
        for comp_gt, count in row.items():
            df.loc[base_gt, comp_gt] = count
    df.index.name = "base_GT"
    df.columns.name = "comp_GT"
    return df


def make_report(summary, label):
    lines = []
    title = "Truvari genotype comparison report" + (f" - {label}" if label else "")
    lines.append(title)
    lines.append("=" * len(title))
    lines.append("")

    lines.append("Site-level matching")
    lines.append("-" * 25)
    lines.append(f"  Base calls total:       {fmt_num(summary.get('base cnt'))}")
    lines.append(f"  Comp calls total:       {fmt_num(summary.get('comp cnt'))}")
    lines.append(f"  TP-base (matched):      {fmt_num(summary.get('TP-base'))}")
    lines.append(f"  TP-comp (matched):      {fmt_num(summary.get('TP-comp'))}")
    lines.append(f"  FP (comp, unmatched):   {fmt_num(summary.get('FP'))}")
    lines.append(f"  FN (base, unmatched):   {fmt_num(summary.get('FN'))}")
    lines.append(f"  Precision:              {fmt_pct(summary.get('precision'))}")
    lines.append(f"  Recall:                 {fmt_pct(summary.get('recall'))}")
    lines.append(f"  F1:                     {fmt_pct(summary.get('f1'))}")
    lines.append("")

    lines.append("Genotype concordance (among site-matched TP calls)")
    lines.append("-" * 52)
    tp_base_tp_gt = summary.get("TP-base_TP-gt")
    tp_base_fp_gt = summary.get("TP-base_FP-gt")
    tp_comp_tp_gt = summary.get("TP-comp_TP-gt")
    tp_comp_fp_gt = summary.get("TP-comp_FP-gt")
    gt_conc = summary.get("gt_concordance")

    lines.append(f"  TP-base with GT match:    {fmt_num(tp_base_tp_gt)}")
    lines.append(f"  TP-base with GT mismatch: {fmt_num(tp_base_fp_gt)}")
    lines.append(f"  TP-comp with GT match:    {fmt_num(tp_comp_tp_gt)}")
    lines.append(f"  TP-comp with GT mismatch: {fmt_num(tp_comp_fp_gt)}")
    lines.append(f"  Overall GT concordance:   {fmt_pct(gt_conc)}")

    # sanity cross-check, if the pieces are all present
    if tp_comp_tp_gt is not None and tp_comp_fp_gt is not None and (tp_comp_tp_gt + tp_comp_fp_gt) > 0:
        derived = tp_comp_tp_gt / (tp_comp_tp_gt + tp_comp_fp_gt)
        lines.append(f"  (derived from TP-comp_TP-gt / (TP-comp_TP-gt + TP-comp_FP-gt): {fmt_pct(derived)})")
    lines.append("")

    return "\n".join(lines)


def make_matrix_section(matrix_df):
    if matrix_df is None:
        return "No gt_matrix present in summary.json.\n"

    lines = ["Genotype confusion matrix (rows = base GT, cols = comp GT)", "-" * 60]
    lines.append(matrix_df.to_string())
    lines.append("")

    row_totals = matrix_df.sum(axis=1)
    lines.append("Per-base-genotype concordance (diagonal / row total):")
    for gt in matrix_df.index:
        total = row_totals[gt]
        if total == 0:
            continue
        match = matrix_df.loc[gt, gt] if gt in matrix_df.columns else 0
        lines.append(f"  base {gt:12s} n={int(total):6d}  concordant={fmt_pct(match / total)}")
    lines.append("")
    return "\n".join(lines)


def make_heatmap(matrix_df, label, outdir):
    if matrix_df is None or matrix_df.empty:
        return

    counts = matrix_df.to_numpy()
    row_totals = counts.sum(axis=1, keepdims=True)
    row_totals_safe = row_totals.copy()
    row_totals_safe[row_totals_safe == 0] = 1
    pct = counts / row_totals_safe * 100

    annot = [[f"{counts[i, j]:,}\n({pct[i, j]:.1f}%)" for j in range(counts.shape[1])]
             for i in range(counts.shape[0])]

    fig, ax = plt.subplots(figsize=(max(5, counts.shape[1] * 1.6), max(4, counts.shape[0] * 1.4)))
    sns.heatmap(pct, annot=annot, fmt="", cmap="Blues", cbar_kws={"label": "% of base GT (row)"},
                xticklabels=matrix_df.columns, yticklabels=matrix_df.index, ax=ax, vmin=0, vmax=100,
                linewidths=0.5, linecolor="white")
    ax.set_xlabel("Comparison genotype")
    ax.set_ylabel("Baseline genotype")
    title = "Genotype confusion matrix" + (f" - {label}" if label else "")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "gt_matrix_heatmap.png"), dpi=150)
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    summary, path = load_summary(args)
    print(f"Loaded {path}")

    matrix_df = build_gt_matrix_df(summary.get("gt_matrix"))

    report = make_report(summary, args.label)
    report += make_matrix_section(matrix_df)

    with open(os.path.join(args.outdir, "gt_stats_report.txt"), "w") as f:
        f.write(report)
    print(report)

    if matrix_df is not None:
        matrix_df.to_csv(os.path.join(args.outdir, "gt_matrix.csv"))
        make_heatmap(matrix_df, args.label, args.outdir)

    print(f"Outputs written to: {args.outdir}")


if __name__ == "__main__":
    main()