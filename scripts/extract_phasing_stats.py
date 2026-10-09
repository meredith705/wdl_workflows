#!/usr/bin/env python3
"""
Extract haplotagging/phasing stats from a margin log file(s) into a CSV.

Usage: python extract_phasing_stats.py margin.log > output.csv
"""

import re
import sys
import csv

def parse_log(text):
    records = []
    current = {}

    # Pattern 1: Separated reads with divisions: H1 X, H2 Y, and H0 Z
    div_re = re.compile(
        r"Separated reads with divisions:\s*H1\s+(\d+),\s*H2\s+(\d+),\s*and\s*H0\s+(\d+)"
    )

    # Pattern 2: Wrote haplotagged BAM in Xm Ys
    bam_time_re = re.compile(
        r"Wrote haplotagged BAM in\s+(\d+)m\s+(\d+)s"
    )

    # Pattern 3: phased VCF / phaseset paths (captures sample-ish identifier)
    vcf_re = re.compile(
        r"Writing phased VCF to (\S+), phaseset info to (\S+)"
    )

    # Pattern 4: variant counts line
    variants_re = re.compile(
        r"Of\s+(\d+)\s+variants:\s+wrote\s+(\d+)\s+with\s+(\d+)\s+phased;\s+"
        r"skipped\s+(\d+)\s+for region,\s+(\d+)\s+for not being analyzed,\s+"
        r"(\d+)\s+for being homozygous,\s+(\d+)\s+for disagreement with margin"
    )

    # Pattern 5: phase sets summary
    phasesets_re = re.compile(
        r"Identified\s+(\d+)\s+phase sets with lengths avg:(-?\d+(?:\.\d+)?),\s*"
        r"min:(-?\d+),\s*max:(-?\d+),\s*N50:(-?\d+)"
    )

    for line in text.splitlines():
        m = div_re.search(line)
        if m:
            if current:
                records.append(current)
            current = {}
            current["H1_reads"] = int(m.group(1))
            current["H2_reads"] = int(m.group(2))
            current["H0_reads"] = int(m.group(3))
            continue

        m = bam_time_re.search(line)
        if m:
            current["bam_write_minutes"] = int(m.group(1))
            current["bam_write_seconds"] = int(m.group(2))
            current["bam_write_total_seconds"] = int(m.group(1)) * 60 + int(m.group(2))
            continue

        m = vcf_re.search(line)
        if m:
            current["phased_vcf_path"] = m.group(1)
            current["phaseset_bed_path"] = m.group(2)
            continue

        m = variants_re.search(line)
        if m:
            current["variants_total"] = int(m.group(1))
            current["variants_wrote"] = int(m.group(2))
            current["variants_phased"] = int(m.group(3))
            current["skipped_region"] = int(m.group(4))
            current["skipped_not_analyzed"] = int(m.group(5))
            current["skipped_homozygous"] = int(m.group(6))
            current["skipped_disagreement"] = int(m.group(7))
            continue

        m = phasesets_re.search(line)
        if m:
            current["phase_sets_count"] = int(m.group(1))
            current["phase_set_len_avg"] = float(m.group(2))
            current["phase_set_len_min"] = int(m.group(3))
            current["phase_set_len_max"] = int(m.group(4))
            current["phase_set_N50"] = int(m.group(5))
            continue

    if current:
        records.append(current)

    return records


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python extract_phasing_stats.py margin.log [margin.log ...]")

    all_records = []
    for fname in sys.argv[1:]:
        with open(fname) as f:
            text = f.read()
        recs = parse_log(text)
        for r in recs:
            r["source_file"] = fname
        all_records.extend(recs)

    if not all_records:
        sys.exit("No matching records found.")

    # Collect all possible fieldnames across records
    fieldnames = []
    for r in all_records:
        for k in r:
            if k not in fieldnames:
                fieldnames.append(k)

    writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames)
    writer.writeheader()
    for r in all_records:
        writer.writerow(r)


if __name__ == "__main__":
    main()