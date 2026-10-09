#!/usr/bin/env python3
"""
quast_summary.py - parse and interpret QUAST reports (report.tsv or report.txt).

Handles multi-assembly reports (one column per assembly, e.g. hap1/hap2),
missing reference metrics ("-"), and prints a side-by-side table plus
plain-language findings tuned for human genome assemblies vs CHM13.

Usage:
    python quast_summary.py quast_out/report.tsv
    python quast_summary.py quast_hap1/ quast_hap2/ report.txt
    python quast_summary.py report.tsv --names hap1,hap2 --tsv summary.tsv
"""
import argparse
import csv
import re
import sys
from pathlib import Path


# ----------------------------------------------------------------- helpers
def num(s):
    """Leading number of a QUAST value ('13 + 177 part' -> 13); None for '-'/blank."""
    if s is None:
        return None
    m = re.match(r"^\s*-?\d+(?:\.\d+)?", s)
    return float(m.group()) if m else None


def fmt_bp(x, signed=False):
    if x is None:
        return "n/a"
    sign = "-" if x < 0 else ("+" if signed and x > 0 else "")
    ax = abs(x)
    for unit, div in (("Gb", 1e9), ("Mb", 1e6), ("kb", 1e3)):
        if ax >= div:
            return f"{sign}{ax / div:.2f} {unit}"
    return f"{sign}{ax:.0f} bp"


def find_report(p):
    p = Path(p)
    if p.is_dir():
        for name in ("report.tsv", "report.txt"):
            if (p / name).exists():
                return p / name
        sys.exit(f"No report.tsv/report.txt in {p}")
    return p


def parse_report(path, fallback_name):
    """Return [(assembly_name, {metric: raw_string}), ...]."""
    header, metrics = None, {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        if "\t" in line:
            parts = line.split("\t")
        else:
            parts = re.split(r"\s{2,}", line.strip())
        parts = [p.strip() for p in parts]
        if parts[0] == "Assembly":
            header = parts[1:]
        elif len(parts) >= 2:
            metrics.setdefault(parts[0], parts[1:])
    n = max((len(v) for v in metrics.values()), default=0)
    if not header or len(header) != n:
        header = [fallback_name] if n == 1 else [f"{fallback_name}:asm{i + 1}" for i in range(n)]
    return [
        (name, {k: (v[i] if i < len(v) else None) for k, v in metrics.items()})
        for i, name in enumerate(header)
    ]


# ---------------------------------------------------------------- analysis
def analyze(d):
    m, findings = {}, []

    def add(level, msg):
        findings.append((level, msg))

    G = lambda k: num(d.get(k))
    total, ref = G("Total length"), G("Reference length")
    aligned, gf = G("Total aligned length"), G("Genome fraction (%)")
    dup = G("Duplication ratio")
    n50, ng50, nga50, n90 = G("N50"), G("NG50"), G("NGA50"), G("N90")
    l90, lg90 = G("L90"), G("LG90")
    has_aln = gf is not None and aligned is not None

    m.update(total=total, ref=ref, dup=dup, gf=gf, n50=n50, ng50=ng50, nga50=nga50,
             contigs=G("# contigs"), mis=G("# misassemblies"),
             mm=G("# mismatches per 100 kbp"), indel=G("# indels per 100 kbp"))

    # --- overall size ------------------------------------------------------
    if total and ref:
        r = total / ref
        m["ratio"] = r
        msg = f"Assembly is {fmt_bp(total)} = {r:.2f}x the reference ({fmt_bp(ref)})."
        if r > 1.8:
            add("WARN", msg + " Close to 2x: looks like both haplotypes concatenated.")
        elif r > 1.25:
            add("WARN", msg + " Much larger than haploid; expect uncollapsed haplotypes/duplication.")
        elif r > 1.08:
            add("NOTE", msg + " Somewhat large; check duplication and unaligned sequence.")
        elif r >= 0.95:
            add("OK", msg + " Close to haploid size.")
        elif r >= 0.85:
            add("NOTE", msg + " Shorter than reference: typical for a haplotype lacking chrY or with "
                "collapsed/missing satellites and rDNA.")
        else:
            add("WARN", msg + " Substantially short; large parts of the genome are likely missing.")

    # --- alignment-based metrics ------------------------------------------
    if not has_aln:
        add("WARN", "No reference-based metrics (all '-'): alignment failed or no reference was given. "
            "Check quast.log / contigs_reports/ and rerun (use --large for human genomes).")
    else:
        # duplication ratio
        if dup is not None:
            if dup <= 1.05:
                add("OK", f"Duplication ratio {dup:.3f}: essentially no redundant sequence.")
            elif dup <= 1.15:
                add("NOTE", f"Duplication ratio {dup:.3f}: minor redundancy.")
            elif dup <= 1.5:
                add("WARN", f"Duplication ratio {dup:.3f}: substantial redundancy (haplotigs/duplicated "
                    "contigs). Consider purge_dups / purge_haplotigs or a higher assembler purge level.")
            else:
                add("WARN", f"Duplication ratio {dup:.3f}: severe; most of the genome is covered more "
                    "than once. Probably unpurged or both haplotypes in one file.")

        # length accounting
        covered = gf / 100.0 * ref
        redundant = max(aligned - covered, 0)
        unal = max(total - aligned, 0)
        uncovered = ref - covered
        m.update(redundant=redundant, unal=unal, uncovered=uncovered)
        excess = total - ref
        msg = (f"Length accounting: {fmt_bp(excess, True)} vs reference = {fmt_bp(redundant)} redundant "
               f"aligned + {fmt_bp(unal)} unaligned - {fmt_bp(uncovered)} reference not covered.")
        if excess > 0 and redundant > 0:
            msg += f" Redundancy explains {100 * redundant / excess:.0f}% of the excess."
        msg += (" (Derived from aligned length minus covered reference; QUAST's duplication ratio is "
                "computed differently, so treat as approximate, and 0 bp means no detectable excess.)")
        add("NOTE", msg)

        # genome fraction
        if gf >= 95:
            add("OK", f"Genome fraction {gf:.2f}%: nearly all of the reference is represented.")
        elif gf >= 90:
            add("NOTE", f"Genome fraction {gf:.2f}% ({fmt_bp(uncovered)} uncovered): likely satellites, "
                "acrocentric short arms/rDNA, or segdups. Intersect uncovered bases with CHM13 censat.")
        else:
            add("WARN", f"Genome fraction {gf:.2f}% ({fmt_bp(uncovered)} uncovered): large missing content.")

        # unaligned
        uf = unal / total if total else 0
        mu = re.match(r"\s*(\d+)(?:\s*\+\s*(\d+)\s*part)?", d.get("# unaligned contigs") or "")
        whole, part = (int(mu.group(1)), int(mu.group(2) or 0)) if mu else (0, 0)
        lvl = "WARN" if uf > 0.05 else ("NOTE" if uf > 0.02 else "OK")
        add(lvl, f"Unaligned: {fmt_bp(unal)} ({100 * uf:.1f}% of assembly); {whole} whole contigs + {part} "
            "partially unaligned. Check *.unaligned.info; screen whole contigs for contamination "
            "(Kraken2/BLAST) and mito/EBV.")

        # misassemblies
        mis, mis_ctg, mis_len = G("# misassemblies"), G("# misassembled contigs"), G("Misassembled contigs length")
        if mis is not None and mis_len is not None and total:
            add("NOTE", f"{mis:.0f} misassemblies in {mis_ctg:.0f} contigs holding {100 * mis_len / total:.0f}% "
                "of the assembly. Against a different individual's reference most are SVs, satellite-length "
                "differences, or duplication artifacts; inspect locations (Icarus / mis_contigs.info) first.")

        # accuracy proxies
        mm, ind = m["mm"], m["indel"]
        if mm is not None:
            if mm > 300:
                add("WARN", f"{mm:.0f} mismatches/100kb: high; check polishing/consensus or ancestry divergence.")
            else:
                add("OK", f"{mm:.0f} mismatches/100kb: in the range expected between two human genomes "
                    "(includes real variation, not just error).")
        if ind is not None and ind > 50:
            add("NOTE", f"{ind:.0f} indels/100kb: on the high side; possible homopolymer errors.")

        # NGA50 vs NG50
        if nga50 and ng50:
            q = nga50 / ng50
            lvl = "OK" if q >= 0.8 else ("NOTE" if q >= 0.5 else "WARN")
            add(lvl, f"NGA50/NG50 = {q:.2f} ({fmt_bp(nga50)} vs {fmt_bp(ng50)}): "
                + ("few breaks at misassemblies." if q >= 0.8 else "contigs get cut at misassemblies."))

    # --- contiguity --------------------------------------------------------
    if ng50 is not None:
        if ng50 >= 30e6:
            tier = "chromosome-arm scale"
        elif ng50 >= 5e6:
            tier = "good"
        elif ng50 >= 1e6:
            tier = "moderate"
        else:
            tier = "fragmented"
        add("OK" if ng50 >= 5e6 else "NOTE",
            f"Contiguity: N50 {fmt_bp(n50)}, NG50 {fmt_bp(ng50)} ({tier} for a human assembly), "
            f"{m['contigs']:.0f} contigs, largest {fmt_bp(G('Largest contig'))}.")
    if n50 and ng50 and total and ref and total > ref and n50 / ng50 < 0.9:
        add("NOTE", "N50 < NG50: excess assembly length is lowering N50, so NG50 is the fairer number.")
    if n90 is not None and n90 < 1e6 and l90 and lg90 and l90 > 3 * lg90:
        add("NOTE", f"Long tail of short contigs (N90 {fmt_bp(n90)}, L90 {l90:.0f} vs LG90 {lg90:.0f}): "
            "often haplotigs/redundant fragments; check their depth and alignment.")

    # --- GC / gaps -----------------------------------------------------------
    gc, rgc = G("GC (%)"), G("Reference GC (%)")
    if gc is not None and rgc is not None and abs(gc - rgc) > 0.4:
        add("NOTE", f"GC {gc:.2f}% vs reference {rgc:.2f}%: noticeable shift; compare content missing "
            "(GC-rich satellites) or extra (AT-rich/contaminant).")
    npk = G("# N's per 100 kbp")
    if npk:
        add("NOTE", f"{npk:.2f} N's per 100 kbp: assembly contains gaps.")
    return m, findings


# ------------------------------------------------------------------ output
TABLE = [
    ("Total length", "total", lambda v: fmt_bp(v)),
    ("x reference", "ratio", lambda v: f"{v:.2f}x"),
    ("Duplication ratio", "dup", lambda v: f"{v:.3f}"),
    ("Genome fraction", "gf", lambda v: f"{v:.2f}%"),
    ("Redundant aligned", "redundant", lambda v: fmt_bp(v)),
    ("Unaligned", "unal", lambda v: fmt_bp(v)),
    ("Reference uncovered", "uncovered", lambda v: fmt_bp(v)),
    ("# contigs", "contigs", lambda v: f"{v:.0f}"),
    ("N50", "n50", lambda v: fmt_bp(v)),
    ("NG50", "ng50", lambda v: fmt_bp(v)),
    ("NGA50", "nga50", lambda v: fmt_bp(v)),
    ("# misassemblies", "mis", lambda v: f"{v:.0f}"),
    ("Mismatches/100kb", "mm", lambda v: f"{v:.1f}"),
    ("Indels/100kb", "indel", lambda v: f"{v:.1f}"),
]
TAG = {"OK": "[ OK ]", "NOTE": "[NOTE]", "WARN": "[WARN]"}


def print_table(results):
    labels = [r[0] for r in results]
    w0 = max(len(t[0]) for t in TABLE)
    w = max(12, max(len(l) for l in labels) + 1)
    print(f"{'':<{w0}}  " + "".join(f"{l:>{w}}" for l in labels))
    for label, key, f in TABLE:
        cells = []
        for _, m, _ in results:
            v = m.get(key)
            cells.append("-" if v is None else f(v))
        print(f"{label:<{w0}}  " + "".join(f"{c:>{w}}" for c in cells))


def main():
    ap = argparse.ArgumentParser(description="Parse and interpret QUAST reports.")
    ap.add_argument("inputs", nargs="+", help="report.tsv / report.txt files or QUAST output dirs")
    ap.add_argument("--names", help="comma-separated assembly names, in column order across all inputs")
    ap.add_argument("--tsv", help="write key metrics for all assemblies to this TSV")
    args = ap.parse_args()

    results = []
    for inp in args.inputs:
        path = find_report(inp)
        fallback = path.parent.name or path.stem
        for name, d in parse_report(path, fallback):
            m, findings = analyze(d)
            results.append((name, m, findings))
    if args.names:
        for (name, m, f), new in zip(list(results), args.names.split(",")):
            results[results.index((name, m, f))] = (new, m, f)

    print_table(results)
    for name, _, findings in results:
        print(f"\n=== {name} ===")
        for level, msg in findings:
            print(f"{TAG[level]} {msg}")

    if len(results) == 2 and all(r[1].get("total") for r in results):
        combined = sum(r[1]["total"] for r in results)
        print(f"\nCombined length of the {len(results)} assemblies: {fmt_bp(combined)} "
              "(a phased diploid human pair should be roughly 5.9-6.2 Gb).")

    if args.tsv:
        keys = [t[1] for t in TABLE]
        with open(args.tsv, "w", newline="") as fh:
            wr = csv.writer(fh, delimiter="\t")
            wr.writerow(["assembly"] + keys)
            for name, m, _ in results:
                wr.writerow([name] + [m.get(k, "") if m.get(k) is not None else "" for k in keys])
        print(f"\nWrote {args.tsv}")


if __name__ == "__main__":
    main()