import pysam, pandas as pd
from collections import defaultdict
import sys, pysam

path = sys.argv[1] if len(sys.argv) > 1 else "-"
vcf = pysam.VariantFile(path)

#vcf = pysam.VariantFile("merged.vcf.gz")

rows = defaultdict(lambda: defaultdict(int))

def size_bin(l):
    l = abs(l)
    return "<100" if l < 100 else "100-1k" if l < 1000 else "1k-10k" if l < 10000 else ">10k"

for rec in vcf:
    svtype = rec.info.get("SVTYPE", "NA")
    svlen = rec.info.get("SVLEN", 0)
    svlen = svlen[0] if isinstance(svlen, tuple) else svlen
    ac = rec.info.get("AC", (0,))[0]
    rare = "singleton" if ac == 1 else "rare" if ac <= 5 else "common"
    for sample, d in rec.samples.items():
        gt = d.get("GT")
        if gt is None or any(a is None for a in gt) or not any(a > 0 for a in gt):
            continue
        r = rows[sample]
        r["total"] += 1
        r[f"type_{svtype}"] += 1
        r[f"size_{size_bin(svlen)}"] += 1
        r[f"af_{rare}"] += 1
        r["hom" if all(a > 0 for a in gt) else "het"] += 1

df = pd.DataFrame(rows).T.fillna(0)
df["het_hom"] = df.het / df.hom.clip(lower=1)

# robust outlier call (median / MAD)
med = df.total.median()
mad = (df.total - med).abs().median()
df["z"] = (df.total - med) / (1.4826 * mad)
df.sort_values("total", ascending=False).to_csv("per_sample_sv_summary.tsv", sep="\t")
print(df[df.z > 3].sort_values("z", ascending=False).head(20))