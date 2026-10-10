#!/usr/bin/env python3
"""
Plot the annotations GATK uses for variant filtering, split into SNPs and indels.

Site-level (INFO/QUAL):  QUAL, QD, FS, SOR, MQ, MQRankSum, ReadPosRankSum, DP
Genotype-level (FORMAT): GQ, DP, and allele balance (from AD, het calls only)

GATK's documented hard-filter thresholds are drawn as dashed lines, and the
fraction of variants failing each one is written to a summary table.
Treat the thresholds as starting points and adjust them for your data.

Usage:
    python plot_gatk_filter_fields.py input.vcf.gz -o gatk_fields
    python plot_gatk_filter_fields.py input.vcf.gz -o gatk_fields --thin 10
    python plot_gatk_filter_fields.py input.vcf.gz -o gatk_fields --region chr20

Outputs:
    <prefix>.png            multi-panel density plots
    <prefix>.summary.tsv    n, median, threshold and % failing per field

Requires: cyvcf2, numpy, pandas, matplotlib
"""
import argparse
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cyvcf2 import VCF

# INFO fields (QUAL is handled separately since it is a column, not INFO)
SITE_FIELDS = ["QUAL", "QD", "FS", "SOR", "MQ", "MQRankSum", "ReadPosRankSum", "DP"]
GT_FIELDS = ["GQ", "DP_gt", "AB"]

# GATK documented hard-filter thresholds: (operator, value) means "fails if x <op> value"
THRESHOLDS = {
    "SNP": {
        "QD": ("<", 2.0), "QUAL": ("<", 30.0), "SOR": (">", 3.0), "FS": (">", 60.0),
        "MQ": ("<", 40.0), "MQRankSum": ("<", -12.5), "ReadPosRankSum": ("<", -8.0),
    },
    "INDEL": {
        "QD": ("<", 2.0), "QUAL": ("<", 30.0), "FS": (">", 200.0),
        "ReadPosRankSum": ("<", -20.0),
    },
}
COLORS = {"SNP": "#3182bd", "INDEL": "#e6550d"}
TITLES = {
    "QUAL": "QUAL", "QD": "QD (QualByDepth)", "FS": "FS (FisherStrand)",
    "SOR": "SOR (StrandOddsRatio)", "MQ": "MQ (RMSMappingQuality)",
    "MQRankSum": "MQRankSum", "ReadPosRankSum": "ReadPosRankSum",
    "DP": "Site DP", "GQ": "Genotype GQ", "DP_gt": "Genotype DP",
    "AB": "Allele balance (het, ALT/(REF+ALT))",
}


def as_float(x):
    if isinstance(x, (tuple, list)):
        x = x[0] if x else None
    try:
        x = float(x)
    except (TypeError, ValueError):
        return None
    return x if np.isfinite(x) else None


def pick(values, mask, rng, k):
    """Valid (>=0) values for masked genotypes, randomly thinned to at most k."""
    vals = values[mask & (values >= 0)]
    if len(vals) > k:
        vals = rng.choice(vals, k, replace=False)
    return vals


def collect(args):
    vcf = VCF(args.vcf, threads=args.threads)
    it = vcf(args.region) if args.region else vcf
    rng = np.random.default_rng(1)
    site = {t: defaultdict(list) for t in COLORS}
    gt = {t: defaultdict(list) for t in COLORS}
    n_seen = 0

    for i, v in enumerate(it):
        if args.thin > 1 and i % args.thin:
            continue
        if v.is_snp:
            vt = "SNP"
        elif v.is_indel:
            vt = "INDEL"
        else:
            continue
        n_seen += 1

        q = as_float(v.QUAL)
        if q is not None:
            site[vt]["QUAL"].append(q)
        for f in SITE_FIELDS[1:]:
            x = as_float(v.INFO.get(f))
            if x is not None:
                site[vt][f].append(x)

        called = v.gt_types != 2  # 2 = missing
        if not called.any():
            continue
        k = args.gt_per_site

        gq = v.format("GQ")
        if gq is not None:
            gt[vt]["GQ"].extend(pick(gq[:, 0].astype(float), called, rng, k))
        dp = v.format("DP")
        if dp is not None:
            gt[vt]["DP_gt"].extend(pick(dp[:, 0].astype(float), called, rng, k))

        # allele balance for het calls at biallelic sites
        if len(v.ALT) == 1:
            ad = v.format("AD")
            if ad is not None and ad.shape[1] >= 2:
                ref, alt = ad[:, 0].astype(float), ad[:, 1].astype(float)
                tot = ref + alt
                ok = (v.gt_types == 1) & (ref >= 0) & (alt >= 0) & (tot > 0)
                if ok.any():
                    ab = np.where(ok, alt / np.where(tot > 0, tot, 1), -1.0)
                    gt[vt]["AB"].extend(pick(ab, ok, rng, k))

    return site, gt, n_seen


def summarize(site, gt, args):
    rows = []
    for vt in COLORS:
        merged = {**site[vt], **{k: v for k, v in gt[vt].items()}}
        for field, vals in merged.items():
            arr = np.asarray(vals, dtype=float)
            if arr.size == 0:
                continue
            thr = THRESHOLDS[vt].get(field)
            if thr is None and field == "GQ":
                thr = ("<", args.gq_min)
            elif thr is None and field == "DP_gt":
                thr = ("<", args.dp_min)
            if thr:
                op, val = thr
                pct = 100 * np.mean(arr < val if op == "<" else arr > val)
            else:
                op, val, pct = "", np.nan, np.nan
            rows.append({"type": vt, "field": field, "n": arr.size,
                         "median": np.median(arr), "fail_op": op,
                         "threshold": val, "pct_failing": pct})
    return pd.DataFrame(rows)


def plot(site, gt, args, n_seen):
    panels = [f for f in SITE_FIELDS + GT_FIELDS
              if any(len(site[t].get(f, [])) or len(gt[t].get(f, [])) for t in COLORS)]
    if not panels:
        raise SystemExit("None of the expected fields were found in the VCF.")

    ncols = 3
    nrows = int(np.ceil(len(panels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.6 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, field in zip(axes, panels):
        data = {t: np.asarray((site[t].get(field) or gt[t].get(field) or []), float)
                for t in COLORS}
        allv = np.concatenate([d for d in data.values() if d.size])
        if field == "AB":
            lo, hi = 0.0, 1.0
        elif field == "GQ":
            lo, hi = 0.0, 99.0
        else:
            lo, hi = np.percentile(allv, [0.5, 99.5])
            if lo == hi:
                hi = lo + 1
        bins = np.linspace(lo, hi, 60)

        for vt, arr in data.items():
            if not arr.size:
                continue
            ax.hist(np.clip(arr, lo, hi), bins=bins, density=True, alpha=0.4,
                    color=COLORS[vt], label=f"{vt} (n={arr.size:,})")
            thr = THRESHOLDS[vt].get(field)
            if thr is None and field == "GQ":
                thr = ("<", args.gq_min)
            if thr is None and field == "DP_gt":
                thr = ("<", args.dp_min)
            if thr and lo <= thr[1] <= hi:
                ax.axvline(thr[1], color=COLORS[vt], ls="--", lw=1.5)
        ax.set_title(TITLES.get(field, field), fontsize=11)
        ax.set_ylabel("Density")
        ax.legend(fontsize=8)

    for ax in axes[len(panels):]:
        ax.axis("off")

    fig.suptitle(f"GATK filtering annotations ({n_seen:,} variants; "
                 f"dashed = threshold; x clipped to 0.5-99.5 percentile)", fontsize=12)
    fig.tight_layout()
    fig.savefig(f"{args.out}.png", dpi=200)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vcf", help="Input VCF/BCF (indexed if using --region)")
    p.add_argument("-o", "--out", default="gatk_fields", help="Output prefix")
    p.add_argument("--region", help="Only read this region, e.g. chr20 or chr20:1-5000000")
    p.add_argument("--thin", type=int, default=1,
                   help="Use every Nth variant (speeds up huge VCFs)")
    p.add_argument("--gt-per-site", type=int, default=5,
                   help="Max genotypes sampled per site for GQ/DP/AB (default 5)")
    p.add_argument("--gq-min", type=float, default=20, help="GQ guide line (default 20)")
    p.add_argument("--dp-min", type=float, default=10, help="Genotype DP guide line (default 10)")
    p.add_argument("--threads", type=int, default=2)
    args = p.parse_args()

    site, gt, n_seen = collect(args)
    if n_seen == 0:
        raise SystemExit("No SNP/indel records found.")

    summary = summarize(site, gt, args)
    summary.to_csv(f"{args.out}.summary.tsv", sep="\t", index=False, float_format="%.4g")
    plot(site, gt, args, n_seen)

    print(summary.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print(f"\nWrote {args.out}.png and {args.out}.summary.tsv")


if __name__ == "__main__":
    main()