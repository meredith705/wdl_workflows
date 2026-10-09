#!/usr/bin/env python3
"""
plot_sv_chroms.py - plot SV locations across chromosomes with genomic annotations.

Each chromosome is a black horizontal line (length = chromosome size) stacked
vertically with chr1 at the bottom and chr22 (or chrX/chrY if requested) at
the top.

SVs are drawn as vertical bars at their positions.

Optional BED annotation regions can be supplied with --regions. The BED file
should contain at least:

    chrom    start    end    type

For example:

    chr1    121500000    124500000    centromere
    chr1    150000000    151000000    segmental_duplication
    chr2    93000000     96000000     centromere

The 4th column determines the annotation type/color.

Centromeres are drawn as thicker lines. Other annotation types are drawn as
thinner lines.

Usage:
    python plot_sv_chroms.py chrom_sizes.bed calls.vcf.gz -o svs.png

    python plot_sv_chroms.py chrom_sizes.bed calls.vcf.gz \
        --regions genomic_regions.bed -o svs_with_regions.png

    python plot_sv_chroms.py chrom_sizes.bed calls.vcf.gz \
        --regions genomic_regions.bed \
        --min-len 8000 \
        --alt-only \
        -o svs_8kb.png
"""

import argparse
import gzip
import re
import sys
from collections import Counter, defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# ---------------------------------------------------------------------------
# SV colors
# ---------------------------------------------------------------------------

COLORS = {
    "DEL": "#d62728",
    "INS": "#1f77b4",
    "DUP": "#2ca02c",
    "INV": "#ff7f0e",
}

OTHER_COLOR = "#7f7f7f"


# ---------------------------------------------------------------------------
# BED annotation colors
#
# Add/change annotation types here as needed.
# Any BED type not listed here gets a color automatically from the fallback
# color list below.
# ---------------------------------------------------------------------------

REGION_COLORS = {
    "centromere": "#9467bd",
    "segmental_duplication": "#8c564b",
    "segdup": "#8c564b",
    "telomere": "#e377c2",
    "repeat": "#bcbd22",
    "gap": "#7f7f7f",
}

# Fallback colors for annotation types not explicitly listed above.
REGION_FALLBACK_COLORS = [
    "#17becf",
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#bcbd22",
    "#7f7f7f",
]


# ---------------------------------------------------------------------------
# File handling
# ---------------------------------------------------------------------------

def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_sizes(path):
    """
    Read chromosome sizes from:
      - .fai
      - 2-column file: chrom length
      - BED: chrom start end
    """
    sizes = {}
    is_fai = str(path).endswith(".fai")

    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue

            f = line.split()

            try:
                if is_fai or len(f) == 2:
                    sizes[f[0]] = int(f[1])
                else:
                    sizes[f[0]] = int(f[2])  # BED end

            except (ValueError, IndexError):
                continue

    return sizes


def read_regions(path):
    """
    Read genomic regions from a BED-like file.

    Expected:
        chrom start end type

    The first three columns are required. The 4th column is the annotation
    type. If no 4th column is present, the region is labeled 'region'.

    Returns:
        list of (chrom, start, end, region_type)
    """
    regions = []

    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue

            f = line.rstrip("\n").split("\t")

            # Allow whitespace-delimited BED files too.
            if len(f) < 3:
                f = line.split()

            if len(f) < 3:
                continue

            try:
                chrom = f[0]
                start = int(f[1])
                end = int(f[2])

                if end <= start:
                    continue

                region_type = f[3].strip() if len(f) >= 4 else "region"

            except (ValueError, IndexError):
                continue

            regions.append((chrom, start, end, region_type))

    return regions


# ---------------------------------------------------------------------------
# VCF parsing
# ---------------------------------------------------------------------------

def parse_info(s):
    d = {}

    for kv in s.split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k] = v
        elif kv:
            d[kv] = True

    return d


def sv_fields(pos, ref, alt, info):
    """Return (svtype, svlen) using INFO tags, falling back to allele lengths."""

    svtype = info.get("SVTYPE")
    svlen = None

    if "SVLEN" in info:
        try:
            svlen = abs(int(float(str(info["SVLEN"]).split(",")[0])))
        except ValueError:
            pass

    first_alt = alt.split(",")[0]

    seq_resolved = (
        not first_alt.startswith("<")
        and not re.search(r"[\[\]]", first_alt)
        and first_alt != "*"
    )

    if svlen is None and seq_resolved:
        svlen = abs(len(first_alt) - len(ref))

    if svlen is None and "END" in info:
        try:
            svlen = abs(int(info["END"]) - pos)
        except ValueError:
            pass

    if svtype is None and seq_resolved:
        d = len(first_alt) - len(ref)

        svtype = (
            "INS"
            if d > 0
            else ("DEL" if d < 0 else "SNV")
        )

    return (svtype or "UNK").upper(), svlen


def carries_alt(fmt, sample):
    keys = fmt.split(":")

    if "GT" not in keys:
        return True

    vals = sample.split(":")
    i = keys.index("GT")

    if i >= len(vals):
        return True

    alleles = re.split(r"[/|]", vals[i])

    return any(a not in ("0", ".") for a in alleles)


def read_vcf(path, min_len, max_len, types, alt_only, pass_only):
    out = []

    with open_text(path) as fh:
        for line in fh:

            if line.startswith("#"):
                continue

            f = line.rstrip("\n").split("\t")

            if len(f) < 8:
                continue

            chrom, pos, _, ref, alt, _, filt, info_s = f[:8]

            if pass_only and filt not in ("PASS", "."):
                continue

            if alt_only and len(f) > 9 and not carries_alt(f[8], f[9]):
                continue

            pos = int(pos)

            svtype, svlen = sv_fields(
                pos,
                ref,
                alt,
                parse_info(info_s),
            )

            if types and svtype not in types:
                continue

            if svlen is None:

                if min_len > 0:
                    continue

            else:

                if svlen < min_len or (
                    max_len and svlen > max_len
                ):
                    continue

            out.append((chrom, pos, svtype, svlen))

    return out


# ---------------------------------------------------------------------------
# Chromosome name handling
# ---------------------------------------------------------------------------

def resolve(chrom, index):
    if chrom in index:
        return chrom

    if "chr" + chrom in index:
        return "chr" + chrom

    if chrom.startswith("chr") and chrom[3:] in index:
        return chrom[3:]

    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():

    ap = argparse.ArgumentParser(
        description="Plot SV positions across chromosomes."
    )

    ap.add_argument(
        "sizes",
        help="chromosome sizes: BED, 2-column, or .fai",
    )

    ap.add_argument(
        "vcf",
        help="VCF/VCF.gz of SVs",
    )

    ap.add_argument(
        "-o",
        "--out",
        default="sv_chromosomes.png",
    )

    # New argument ----------------------------------------------------------

    ap.add_argument(
        "--regions",
        help=(
            "optional BED file of genomic annotations. "
            "Columns: chrom, start, end, type"
        ),
    )

    ap.add_argument(
        "--include-xy",
        action="store_true",
        help="also plot chrX and chrY above chr22",
    )

    ap.add_argument(
        "--chroms",
        help="comma-separated chromosome order, bottom to top "
             "(overrides default)",
    )

    ap.add_argument(
        "--min-len",
        type=int,
        default=0,
        help="minimum SV length (bp)",
    )

    ap.add_argument(
        "--max-len",
        type=int,
        default=0,
        help="maximum SV length (bp), 0 = no limit",
    )

    ap.add_argument(
        "--types",
        help="comma-separated SVTYPEs to keep, e.g. INS,DEL",
    )

    ap.add_argument(
        "--alt-only",
        action="store_true",
        help="only records where the first sample carries an alt allele",
    )

    ap.add_argument(
        "--pass-only",
        action="store_true",
        help="only FILTER == PASS or '.'",
    )

    ap.add_argument(
        "--single-color",
        action="store_true",
        help="draw all bars in one color instead of by SVTYPE",
    )

    ap.add_argument(
        "--bar-height",
        type=float,
        default=0.7,
        help="SV bar height in chromosome-spacing units",
    )

    ap.add_argument(
        "--bar-width",
        type=float,
        default=0.6,
        help="SV bar line width in points",
    )

    # New annotation controls ----------------------------------------------

    ap.add_argument(
        "--region-width",
        type=float,
        default=6.0,
        help="width of ordinary BED annotation lines in points",
    )

    ap.add_argument(
        "--centromere-width",
        type=float,
        default=12.0,
        help="width of centromere annotation lines in points",
    )

    ap.add_argument(
        "--region-alpha",
        type=float,
        default=0.8,
        help="alpha/transparency of BED annotations",
    )

    ap.add_argument(
        "--alpha",
        type=float,
        default=0.7,
    )

    ap.add_argument(
        "--width",
        type=float,
        default=12,
    )

    ap.add_argument(
        "--height",
        type=float,
        default=9,
    )

    ap.add_argument(
        "--dpi",
        type=int,
        default=200,
    )

    ap.add_argument(
        "--title",
    )

    args = ap.parse_args()

    # ----------------------------------------------------------------------
    # Read chromosome sizes
    # ----------------------------------------------------------------------

    sizes = read_sizes(args.sizes)

    if args.chroms:

        order = args.chroms.split(",")

    else:

        order = [f"chr{i}" for i in range(1, 23)]

        if args.include_xy:
            order += ["chrX", "chrY"]

    # Allow chr1 vs 1 naming differences between the sizes file
    # and the requested chromosome list.

    chroms = []

    for c in order:

        r = resolve(c, sizes)

        if r is None:

            print(
                f"warning: {c} not in size file, skipping",
                file=sys.stderr,
            )

        else:

            chroms.append((c, r))

    if not chroms:
        sys.exit(
            "None of the requested chromosomes were found in the size file."
        )

    ypos = {
        real: i
        for i, (_, real) in enumerate(chroms)
    }

    # ----------------------------------------------------------------------
    # Read VCF
    # ----------------------------------------------------------------------

    types = (
        {t.strip().upper() for t in args.types.split(",")}
        if args.types
        else None
    )

    calls = read_vcf(
        args.vcf,
        args.min_len,
        args.max_len,
        types,
        args.alt_only,
        args.pass_only,
    )

    by_type = defaultdict(lambda: ([], []))

    counts, skipped = Counter(), 0

    for chrom, pos, svtype, _ in calls:

        real = resolve(chrom, ypos)

        if real is None:

            skipped += 1
            continue

        key = "SV" if args.single_color else svtype

        by_type[key][0].append(pos / 1e6)
        by_type[key][1].append(ypos[real])

        counts[key] += 1

    # ----------------------------------------------------------------------
    # Read optional BED annotations
    # ----------------------------------------------------------------------

    regions = []

    if args.regions:

        regions = read_regions(args.regions)

    # Determine colors for annotation types.

    region_types = []

    for _, _, _, region_type in regions:

        if region_type not in region_types:
            region_types.append(region_type)

    region_color_map = {}

    fallback_index = 0

    for region_type in region_types:

        # Explicit color takes priority.

        if region_type.lower() in REGION_COLORS:

            region_color_map[region_type] = (
                REGION_COLORS[region_type.lower()]
            )

        else:

            region_color_map[region_type] = (
                REGION_FALLBACK_COLORS[
                    fallback_index % len(REGION_FALLBACK_COLORS)
                ]
            )

            fallback_index += 1

    # ----------------------------------------------------------------------
    # Plot
    # ----------------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(args.width, args.height)
    )

    # ----------------------------------------------------------------------
    # Chromosome lines
    # ----------------------------------------------------------------------

    for _, real in chroms:

        ax.hlines(
            ypos[real],
            0,
            sizes[real] / 1e6,
            color="black",
            lw=2.5,
            zorder=3,
        )

    # ----------------------------------------------------------------------
    # BED genomic annotations
    #
    # These are drawn BEFORE the SV bars so the SV calls remain visible.
    # ----------------------------------------------------------------------

    region_handles = []

    if regions:

        for chrom, start, end, region_type in regions:

            real = resolve(chrom, ypos)

            if real is None:
                continue

            color = region_color_map[region_type]

            # Centromeres get substantially thicker lines.
            if region_type.lower() == "centromere":
                linewidth = args.centromere_width
                zorder = 4
            else:
                linewidth = args.region_width
                zorder = 3.5

            ax.plot(
                [start / 1e6, end / 1e6],
                [ypos[real], ypos[real]],
                color=color,
                lw=linewidth,
                alpha=args.region_alpha,
                zorder=zorder,
                solid_capstyle="butt",
            )

        # Separate legend for BED annotations
        for region_type in region_types:

            color = region_color_map[region_type]

            linewidth = (
                args.centromere_width
                if region_type.lower() == "centromere"
                else args.region_width
            )

            region_handles.append(
                Line2D(
                    [0],
                    [0],
                    color=color,
                    lw=linewidth,
                    label=region_type,
                )
            )

    # ----------------------------------------------------------------------
    # SV bars
    # ----------------------------------------------------------------------

    sv_handles = []

    for key in sorted(
        by_type,
        key=lambda k: -counts[k],
    ):

        xs, ys = by_type[key]

        color = (
            OTHER_COLOR
            if args.single_color
            else COLORS.get(key, OTHER_COLOR)
        )

        ax.vlines(
            xs,
            ys,
            [y + args.bar_height for y in ys],
            color=color,
            lw=args.bar_width,
            alpha=args.alpha,
            zorder=5,
        )

        sv_handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                lw=3,
                label=f"{key} (n={counts[key]})",
            )
        )

    # ----------------------------------------------------------------------
    # Axes
    # ----------------------------------------------------------------------

    ax.set_yticks(
        range(len(chroms))
    )

    ax.set_yticklabels(
        [c for c, _ in chroms]
    )

    ax.set_ylim(
        -0.5,
        len(chroms) - 0.1,
    )

    ax.set_xlim(
        0,
        max(
            sizes[r]
            for _, r in chroms
        ) / 1e6 * 1.01,
    )

    ax.set_xlabel(
        "Position (Mb)"
    )

    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    ax.set_title(
        args.title
        or f"8KB or larger SV locations (n={sum(counts.values())})"
    )

    # ----------------------------------------------------------------------
    # Legends
    # ----------------------------------------------------------------------

    # SV legend on the upper right.

    if sv_handles:

        legend_sv = ax.legend(
            handles=sv_handles,
            loc="upper right",
            frameon=False,
            title="SV type",
        )

        ax.add_artist(legend_sv)

    # BED annotation legend on the lower/right side.

    if region_handles:

        ax.legend(
            handles=region_handles,
            loc="center right",
            frameon=False,
            title="Genomic annotation",
        )

    # ----------------------------------------------------------------------
    # Save
    # ----------------------------------------------------------------------

    fig.tight_layout()

    fig.savefig(
        args.out,
        dpi=args.dpi,
    )

    print(
        f"Plotted {sum(counts.values())} variants "
        f"on {len(chroms)} chromosomes -> {args.out}"
    )

    if regions:

        print(
            f"Plotted {len(regions)} genomic annotation regions "
            f"from {args.regions}"
        )

    if skipped:

        print(
            f"Skipped {skipped} variants on chromosomes "
            f"not in the plot.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()