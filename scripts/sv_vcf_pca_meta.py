#!/usr/bin/env python3
"""
PCA of a structural-variant VCF, with samples colored by diagnosis (or cohort).

Cohort is taken from the sample-name prefix, e.g. with the default delimiter
"HGSVC_NA19240" -> "HGSVC", "1KG-HG00096" -> "1KG".

Diagnosis (or any other metadata column) is read from a metadata table
(TSV/CSV, delimiter auto-detected) passed with --meta. The table must have a
column with sample names (--meta-sample-col) matching the VCF sample IDs, and a
column with the diagnosis (--meta-diagnosis-col).

Genotypes are converted to alt-allele dosage (0/1/2), variants are filtered
on missingness and allele frequency, missing values are mean-imputed, and
each variant is standardized before PCA.

Usage:
    python sv_vcf_pca.py -i svs.vcf.gz -m meta.tsv -o sv_pca
    python sv_vcf_pca.py -i svs.vcf.gz -m meta.csv --meta-sample-col IID --meta-diagnosis-col Dx
    python sv_vcf_pca.py -i svs.vcf.gz -m meta.tsv --color-by cohort
    python sv_vcf_pca.py -i svs.vcf.gz -o sv_pca --prefix-len 3   # no metadata -> colored by cohort

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

# Fixed cohort colors (any cohort not listed falls back to a tab20 color)
COHORT_COLORS = {
    "PPMI": "#2ca02c",
    "RUSH": "#d62728",
    "HBCC": "#ff7f0e",
    "NABEC": "#1f77b4",
}
UNKNOWN_LABEL = "Unknown"
UNKNOWN_COLOR = "#b0b0b0"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--vcf", required=True, help="Input VCF/BCF (.vcf, .vcf.gz, .bcf)")
    p.add_argument("-m", "--meta", default=None,
                   help="Metadata table (TSV/CSV, delimiter auto-detected) with sample IDs and diagnosis")
    p.add_argument("--meta-sample-col", default=None,
                   help="Metadata column with sample IDs (default: first column)")
    p.add_argument("--meta-diagnosis-col", default="diagnosis",
                   help="Metadata column with diagnosis (default: 'diagnosis')")
    p.add_argument("--color-by", choices=["diagnosis", "cohort"], default=None,
                   help="Color points by diagnosis or cohort (default: diagnosis if --meta given, else cohort)")
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
    args = p.parse_args()

    if args.color_by is None:
        args.color_by = "diagnosis" if args.meta else "cohort"
    if args.color_by == "diagnosis" and not args.meta:
        p.error("--color-by diagnosis requires --meta")
    return args


def get_cohort(sample, delimiter, prefix_len):
    if prefix_len:
        return sample[:prefix_len]
    return re.split(delimiter, sample, maxsplit=1)[0]


def load_diagnosis(meta_path, samples, sample_col, diag_col):
    """Return a list of diagnosis labels aligned to `samples` (UNKNOWN_LABEL if missing)."""
    meta = pd.read_csv(meta_path, sep=None, engine="python", dtype=str)
    meta.columns = [c.strip() for c in meta.columns]

    if sample_col is None:
        sample_col = meta.columns[0]
        print(f"--meta-sample-col not given; using first column '{sample_col}'", file=sys.stderr)
    for col in (sample_col, diag_col):
        if col not in meta.columns:
            sys.exit(f"Column '{col}' not found in {meta_path}. Available columns: {list(meta.columns)}")

    meta[sample_col] = meta[sample_col].str.strip()
    meta[diag_col] = meta[diag_col].str.strip()
    dups = meta[sample_col].duplicated(keep=False)
    if dups.any():
        print(f"WARNING: {meta.loc[dups, sample_col].nunique()} duplicated sample IDs in metadata; "
              f"keeping the first occurrence of each", file=sys.stderr)
        meta = meta.drop_duplicates(subset=sample_col, keep="first")

    m2 = meta.set_index(sample_col)
    print(" meta cols:\n",meta.columns)
    print(" m2:\n", m2.columns)
    print(m2.head())

    lookup = meta.set_index(sample_col)[diag_col]
    diagnosis = []
    missing = []
    for s in samples:
        val = lookup.get(s)
        if val is None or pd.isna(val) or val == "":
            diagnosis.append(UNKNOWN_LABEL)
            missing.append(s)
        else:
            diagnosis.append(val)

    if missing:
        shown = ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else "")
        print(f"WARNING: {len(missing)} / {len(samples)} VCF samples have no diagnosis in metadata "
              f"(labelled '{UNKNOWN_LABEL}'): {shown}", file=sys.stderr)
    extra = set(meta[sample_col]) - set(samples)
    if extra:
        print(f"Note: {len(extra)} metadata samples are not in the VCF", file=sys.stderr)
    return diagnosis


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


def build_color_map(groups, color_by):
    """Map each group label to a color."""
    cmap = plt.get_cmap("tab10" if len(groups) <= 10 else "tab20")
    colors = {}
    i = 0
    for g in groups:
        if g == UNKNOWN_LABEL:
            colors[g] = UNKNOWN_COLOR
        elif color_by == "cohort" and g in COHORT_COLORS:
            colors[g] = COHORT_COLORS[g]
        else:
            colors[g] = cmap(i % cmap.N)
            i += 1
    return colors


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
    if args.meta:
        df.insert(0, "diagnosis", load_diagnosis(args.meta, samples, args.meta_sample_col,
                                                 args.meta_diagnosis_col))
    df.insert(0, "sample", samples)
    tsv = f"{args.out_prefix}_pcs.tsv"
    df.to_csv(tsv, sep="\t", index=False)
    print(f"Wrote {tsv}", file=sys.stderr)
    print("Cohort sizes:\n" + df["cohort"].value_counts().to_string(), file=sys.stderr)
    if args.meta:
        print("Diagnosis sizes:\n" + df["diagnosis"].value_counts().to_string(), file=sys.stderr)

    # Plot
    px, py = args.pcs
    color_col = args.color_by
    # Put "Unknown" last so it's drawn on top of nothing and listed at the bottom of the legend
    groups = sorted(g for g in df[color_col].unique() if g != UNKNOWN_LABEL)
    if UNKNOWN_LABEL in set(df[color_col]):
        groups.append(UNKNOWN_LABEL)
    color_of = build_color_map(groups, color_col)

    fig, ax = plt.subplots(figsize=(8, 6.5))
    for g in groups:
        sub = df[df[color_col] == g]
        ax.scatter(sub[f"PC{px}"], sub[f"PC{py}"], s=args.point_size, color=color_of[g],
                   label=f"{g} (n={len(sub)})", alpha=0.85, edgecolor="k", linewidth=0.3)
    if args.label:
        for _, r in df.iterrows():
            ax.annotate(r["sample"], (r[f"PC{px}"], r[f"PC{py}"]), fontsize=6, alpha=0.7,
                        xytext=(3, 3), textcoords="offset points")

    ax.set_xlabel(f"PC{px} ({var[px - 1]:.1f}%)")
    ax.set_ylabel(f"PC{py} ({var[py - 1]:.1f}%)")
    ax.set_title(f"{args.plot_title_prefix} SV genotype PCA")
    ax.axhline(0, color="grey", lw=0.5, zorder=0)
    ax.axvline(0, color="grey", lw=0.5, zorder=0)
    ax.legend(title=color_col.capitalize(), bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
    fig.tight_layout()

    out_plot = f"{args.out_prefix}_PC{px}_PC{py}_by_{color_col}.{args.format}"
    fig.savefig(out_plot, dpi=300, bbox_inches="tight")
    print(f"Wrote {out_plot}", file=sys.stderr)


if __name__ == "__main__":
    main()