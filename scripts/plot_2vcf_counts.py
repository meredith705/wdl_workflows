import argparse
import pandas as pd
import pysam
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import datetime

"""
	Script to count up SV types and VCF alleles per sample. 

	Usage: python3 vcf_variant_counts.py -i your.vcf.gz 

	Add: --plot_violin to output a violin+swarm plot of the counts per sample

	Author: Melissa Meredith UCSC
	02/2025
"""

labelsize=16
ticksize=14

sns.set(
    rc={
        "xtick.labelsize": labelsize,
        "ytick.labelsize": labelsize,
        "axes.labelsize": labelsize,
        "axes.titlesize": ticksize
    }
)



if __name__ == "__main__":

	# Make an argument parser
	parser = argparse.ArgumentParser(description="Process a vcf file using pysam.")

	# add argument for the input vcf file
	parser.add_argument(
		"-a","--in_hifi_count_file",
		type=str,
		required=True,
		help="Path to the input vcf file to be analyzed. Can be bgzipped, having an index would increase processing speed."
	)

	parser.add_argument(
		"-b","--in_shasta_count_file",
		type=str,
		required=True,
		help="Path to the input vcf file to be analyzed. Can be bgzipped, having an index would increase processing speed."
	)

	args = parser.parse_args()

	hifi = pd.read_csv(args.in_hifi_count_file, sep="\t")
	hifi['assembler'] = "hifiasm"
	shasta = pd.read_csv(args.in_shasta_count_file, sep="\t")
	shasta['assembler'] = "shasta"

	# print('hifi', hifi.head())
	# print('shasta', shasta.head())

	data = pd.concat([hifi, shasta])

	# plot side by side barchart
	plt.figure(figsize=(10,8))

	ax = sns.barplot(
		data=data,
		x="type",
		y="count",
		hue="assembler"
	)
	# Add labels on bars
	for p in ax.patches:
	    height = p.get_height()
	    ax.text(
	        p.get_x() + p.get_width() / 2,
	        height,
	        f"{int(height)}",
	        ha="center",
	        va="bottom",
	        fontsize=10
	    )

	plt.xlabel("SV type")
	plt.ylabel("Count")
	plt.title("SV counts by assembler")
	plt.xticks(rotation=45)

	plt.tight_layout()

	prefix = args.in_shasta_count_file.split(".")[0].split("_hapdiff")[0]
	plt.savefig(prefix+"vcf_barchart.png")


