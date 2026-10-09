#!/usr/bin/env python3
"""
Rewrite the ID column of an SV VCF as:

    svim_asm_<SVTYPE>_<abs(SVLEN)>_<AC>[_<CollapseId>]

CollapseId is only appended when it is present in INFO.

Usage:
    python rename_sv_ids.py input.vcf.gz | bgzip > renamed.vcf.gz
    python rename_sv_ids.py input.vcf --unique > renamed.vcf

Notes:
  - SVLEN missing (e.g. BND): falls back to abs(END - POS), else "NA".
  - AC missing: "NA". Multi-valued AC (multiallelic) is joined with "-".
  - --unique appends _2, _3, ... to any ID that would otherwise repeat.
Reads plain or gzipped VCF; writes plain VCF to stdout.
"""
import argparse
import gzip
import sys
from collections import defaultdict


def open_text(path):
    if path == "-":
        return sys.stdin
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path)


def parse_info(info):
    d = {}
    if info == ".":
        return d
    for item in info.split(";"):
        k, _, v = item.partition("=")
        d[k] = v  # flags get an empty string
    return d


def first(value):
    return value.split(",")[0]


def abs_len(info, pos):
    if "SVLEN" in info and info["SVLEN"] not in ("", "."):
        try:
            return str(abs(int(float(first(info["SVLEN"])))))
        except ValueError:
            pass
    if "END" in info:
        try:
            return str(abs(int(info["END"]) - int(pos)))
        except ValueError:
            pass
    return "NA"


def make_id(fields):
    info = parse_info(fields[7])
    svtype = info.get("SVTYPE") or "NA"
    length = abs_len(info, fields[1])
    ac = info.get("AC", "")
    ac = ac.replace(",", "-") if ac not in ("", ".") else "NA"

    parts = ["svim_asm", svtype, length, ac]
    cid = info.get("CollapseId", "")
    if cid not in ("", "."):
        parts.append(cid)
    return "_".join(parts)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vcf", help="Input VCF/VCF.gz ('-' for stdin)")
    p.add_argument("--unique", action="store_true",
                   help="Append _N to duplicate IDs so all IDs are unique")
    args = p.parse_args()

    seen = defaultdict(int)
    out = sys.stdout
    for line in open_text(args.vcf):
        if line.startswith("#"):
            out.write(line)
            continue
        fields = line.rstrip("\n").split("\t")
        new_id = make_id(fields)
        if args.unique:
            seen[new_id] += 1
            if seen[new_id] > 1:
                new_id = f"{new_id}_{seen[new_id]}"
        fields[2] = new_id
        out.write("\t".join(fields) + "\n")


if __name__ == "__main__":
    main()