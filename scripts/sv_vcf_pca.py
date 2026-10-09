#!/usr/bin/env python3
"""
PCA of a structural-variant VCF, with samples colored by cohort.

Cohort is taken from the sample-name prefix, e.g. with the default delimiter
"HGSVC_NA19240" -> "HGSVC", "1KG-HG00096" -> "1KG".

Genotypes are converted to alt-allele dosage (0/1/2), variants are filtered
on missingness and allele frequency, missing values are mean-imputed, and
each variant is standardized before PCA.

Usage:
    python sv_vcf_pca.py -i svs.vcf.gz -o sv_pca
    python sv_vcf_pca.py -i svs.vcf.gz -o sv_pca --prefix-len 3
    python sv_vcf_pca.py -i svs.vcf.gz -o sv_pca --delimiter "_" --pcs 2 3

Requires: pysam numpy pandas scikit-learn matplotlib
"""

import argparse
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
from sklearn.decomposition import PCA


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--vcf", required=True, help="Input VCF/BCF (.vcf, .vcf.gz, .bcf)")
    p.add_argument("-o", "--out-prefix", default="sv_pca", help="Output prefix (default: sv_pca)")
    p.add_argument("--delimiter", default=r"[_\-.]",
                   help=r"Regex separating cohort prefix from the rest of the sample name (default: '[_\-.]')")
    p.add_argument("--prefix-len", type=int, default=None,
                   help="Use the first N characters of the sample name as cohort (overrides --delimiter)")
    p.add_argument("--pass-only", action="store_true", help="Keep only PASS (or unfiltered '.') records")
    p.add_argument("--max-missing", type=float, default=0.1,
                   help="Drop variants with missing-genotype rate above this (default: 0.1)")
    p.add_argument("--min-af", type=float, default=0.01,
                   help="Drop variants with alt allele frequency below this or above 1 - this (default: 0.01)")
    p.add_argument("--n-pcs", type=int, default=10, help="Number of PCs to compute (default: 10)")
    p.add_argument("--pcs", type=int, nargs=2, default=[1, 2], metavar=("X", "Y"),
                   help="Which PCs to plot (1-based, default: 1 2)")
    p.add_argument("--no-scale", action="store_true", help="Center only; don't scale variants to unit variance")
    p.add_argument("--label", action="store_true", help="Annotate each point with its sample name")
    p.add_argument("--point-size", type=float, default=40)
    p.add_argument("--format", default="png", choices=["png", "pdf", "svg"], help="Plot format (default: png)")
    p.add_argument("--plot_title_prefix", type=str, default="SV")
    return p.parse_args()


def get_cohort(sample, delimiter, prefix_len):
    if prefix_len:
        return sample[:prefix_len]
    return re.split(delimiter, sample, maxsplit=1)[0]


def load_dosage_matrix(vcf_path, pass_only=False):
    """Return (samples, matrix[n_samples x n_variants]) of alt-allele dosage; NaN = missing."""
    columns = []
    with pysam.VariantFile(vcf_path) as vcf:
        samples = list(vcf.header.samples)
        n_total = 0
        for rec in vcf:
            n_total += 1
            if pass_only and len(rec.filter) and "PASS" not in rec.filter.keys():
                continue
            col = np.full(len(samples), np.nan, dtype=np.float32)
            for i, s in enumerate(samples):
                gt = rec.samples[s]["GT"]
                if gt is None or len(gt) == 0 or None in gt:
                    continue
                col[i] = sum(1 for a in gt if a > 0)
            columns.append(col)

    print(f"Read {n_total} records; {len(columns)} kept after PASS filter", file=sys.stderr)
    if not columns:
        sys.exit("No variants loaded.")
    return samples, np.column_stack(columns)


def filter_and_impute(X, max_missing, min_af):
    missing_rate = np.isnan(X).mean(axis=0)
    with np.errstate(all="ignore"):
        af = np.nanmean(X, axis=0) / 2.0  # assumes diploid
    keep = (missing_rate <= max_missing) & (af >= min_af) & (af <= 1 - min_af)
    print(f"{keep.sum()} / {X.shape[1]} variants pass missingness/AF filters", file=sys.stderr)
    if keep.sum() < 2:
        sys.exit("Too few variants left after filtering; relax --max-missing / --min-af.")
    X = X[:, keep].copy()

    # Mean-impute missing genotypes per variant
    col_means = np.nanmean(X, axis=0)
    nan_r, nan_c = np.where(np.isnan(X))
    X[nan_r, nan_c] = col_means[nan_c]
    return X


def main():
    args = parse_args()

    samples, X = load_dosage_matrix(args.vcf, args.pass_only)
    X = filter_and_impute(X, args.max_missing, args.min_af)

    # Center (and optionally scale) each variant
    X = X - X.mean(axis=0)
    if not args.no_scale:
        sd = X.std(axis=0)
        sd[sd == 0] = 1.0
        X = X / sd

    n_pcs = min(args.n_pcs, X.shape[0], X.shape[1])
    if max(args.pcs) > n_pcs:
        sys.exit(f"Requested PC {max(args.pcs)} but only {n_pcs} were computed.")
    pca = PCA(n_components=n_pcs)
    coords = pca.fit_transform(X)
    var = pca.explained_variance_ratio_ * 100

    cohorts = [get_cohort(s, args.delimiter, args.prefix_len) for s in samples]
    df = pd.DataFrame(coords, columns=[f"PC{i + 1}" for i in range(n_pcs)])
    df.insert(0, "cohort", cohorts)
    df.insert(0, "sample", samples)
    tsv = f"{args.out_prefix}_pcs.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    print(f"Wrote {tsv}", file=sys.stderr)
    print("Cohort sizes:\n" + df["cohort"].value_counts().to_string(), file=sys.stderr)

    # Plot
    px, py = args.pcs
    unique = sorted(df["cohort"].unique())
    cmap = plt.get_cmap("tab10" if len(unique) <= 10 else "tab20")
    # color_of = {c: cmap(i % cmap.N) for i, c in enumerate(unique)}
    color_of = {"PPMI" : '#2ca02c',
                     "RUSH" : '#d62728',
                     "HBCC" : '#ff7f0e',
                     "NABEC" : '#1f77b4'
                    }

    fig, ax = plt.subplots(figsize=(8, 6.5))
    for c in unique:
        sub = df[df["cohort"] == c]
        ax.scatter(sub[f"PC{px}"], sub[f"PC{py}"], s=args.point_size, color=color_of[c],
                   label=f"{c} (n={len(sub)})", alpha=0.85, edgecolor="k", linewidth=0.3)
    if args.label:
        for _, r in df.iterrows():
            ax.annotate(r["sample"], (r[f"PC{px}"], r[f"PC{py}"]), fontsize=6, alpha=0.7,
                        xytext=(3, 3), textcoords="offset points")

    ax.set_xlabel(f"PC{px} ({var[px - 1]:.1f}%)")
    ax.set_ylabel(f"PC{py} ({var[py - 1]:.1f}%)")
    ax.set_title(f"{args.plot_title_prefix} SV genotype PCA")
    ax.axhline(0, color="grey", lw=0.5, zorder=0)
    ax.axvline(0, color="grey", lw=0.5, zorder=0)
    ax.legend(title="Cohort", bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()

    out_plot = f"{args.out_prefix}_PC{px}_PC{py}.{args.format}"
    fig.savefig(out_plot, dpi=300, bbox_inches="tight")
    print(f"Wrote {out_plot}", file=sys.stderr)


if __name__ == "__main__":
    main()