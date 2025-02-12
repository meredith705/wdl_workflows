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

def vcfEntriesPerSample(in_vcf):
	"""
	Function that takes in a vcf with multiple samples and counts up genotypes of vcf entries per sample

	for each sample column the number of alleles is counted and output as a tsv. 

	"""



	# # Open the VCF file using pysam
	vcf_file = pysam.VariantFile(in_vcf) 

	# set up a dictionary to track the number of variant allels per sample 
	variant_counts = {sample:0 for sample in vcf_file.header.samples}
	# make a dictionary of SV types in the vcf
	svTypes = {}
	# keep track of the number of variants in the vcf
	num_records = 0

	for record in vcf_file:

		# increment the record counter
		num_records+=1

		# add the SV type to the dictionary and/or increment the count and keep track of the SV length
		if record.info['SVTYPE'] not in svTypes.keys():
			svTypes[record.info['SVTYPE']]={'count':1,'lenghts':[record.info['SVLEN']]}
		else:
			svTypes[record.info['SVTYPE']]['count'] += 1
			svTypes[record.info['SVTYPE']]['lenghts'].append(record.info['SVLEN'])

		# for each sample entry in the vcf record add up the alleles ( 1's )
		for sample, data in record.samples.items():
			
			if 'GT' in data:
				# print(record.info['SVTYPE'], 'sample', sample, data["GT"])
				for allele in data["GT"]:
					# count alleles that are not '.' nor 0
					if allele is not None and allele > 0:
						# increment the variant count for the sample 
						variant_counts[sample] += 1

	vcf_file.close()


	print(f'finished analyzing VCF: {num_records} variants in the file' )
	print(svTypes)

	# write out the variant counts to a tsv file, including the header column names, but not the index
	variant_count_df = pd.DataFrame( list(variant_counts.items()), columns=['Sample', 'VarintCount'])
	variant_count_df.to_csv("sample_variant_counts.tsv", header=True, index=False, sep="\t")


	# return the variant count dataframe to main
	return variant_count_df

def violin_swarm(x,y,data,ax,swarm_pt_size = 3):
	""" make a violin plot with a swarm of datapoints on top """

	sns.violinplot(x=x, y=y, data=data, cut=0.25, inner="quartile", alpha = 0.05, ax=ax, edgecolor='black')
	sns.swarmplot(x=x, y=y, data=data, s=swarm_pt_size, alpha=1, ax=ax, color='black')


def plot_violin(vcf_data):
	""" set up the figure to plot a violin """
	fig, axs = plt.subplots(figsize=(8,4))

	violin_swarm(['samples']*vcf_data.shape[0], 'VarintCount', vcf_data, axs)

	plt.tight_layout()
	plt.savefig("sample_variant_counts.png", dpi=300)




if __name__ == "__main__":

	# Make an argument parser
	parser = argparse.ArgumentParser(description="Process a vcf file using pysam.")

	# add argument for the input vcf file
	parser.add_argument(
		"-i","--in_vcf_file",
		type=str,
		required=True,
		help="Path to the input vcf file to be analyzed. Can be bgzipped, having an index would increase processing speed."
	)

	# add arugment for making a plot of sample variant counts 
	parser.add_argument(
		'--plot_violin', 
		action='store_true', 
		help='make a violin plot of the number of variants per sample'
	)

	if len(sys.argv) == 0:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	args = parser.parse_args()

	# Process the VCF file
	vcf_df = vcfEntriesPerSample(args.in_vcf_file)

	if args.plot_violin:

		plot_violin(vcf_df)	


