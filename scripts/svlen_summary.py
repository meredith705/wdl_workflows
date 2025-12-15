#!/usr/bin/env python3
import argparse
import gzip

def open_vcf(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def sum_ins_del(vcf_path):
    total_ins = 0
    total_del = 0

    with open_vcf(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            ref = fields[3]
            alts = fields[4].split(",")

            for alt in alts:
                # skip symbolic calls like <INS>, <DEL>, <DUP>, etc
                if alt.startswith("<") and alt.endswith(">"):
                    continue

                # compute length delta
                d = len(alt) - len(ref)

                if d > 0:
                    total_ins += d
                elif d < 0:
                    total_del += abs(d)

    return total_ins, total_del


def main():
    parser = argparse.ArgumentParser(description="Sum INS and DEL sequence length from a VCF")
    parser.add_argument("-v","--vcf", help="Input VCF or VCF.gz")

    args = parser.parse_args()

    total_ins, total_del = sum_ins_del(args.vcf)

    print(f"Total inserted sequence: {total_ins}")
    print(f"Total deleted sequence:  {total_del}")


if __name__ == "__main__":
    main()
