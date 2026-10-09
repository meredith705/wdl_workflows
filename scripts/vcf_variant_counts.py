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

labelsize=14
ticksize=24

sns.set(style="whitegrid")

sns.set(
    rc={
        "xtick.labelsize": labelsize,
        "ytick.labelsize": labelsize,
        "axes.labelsize": labelsize,
        "axes.titlesize": ticksize
    }
)



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
			if 'SVLEN' not in record.info:
				ref_len = len(record.ref)
				alt_len = len(record.alts)
				if alt_len - ref_len == 0:
					varLen = ref_len
				else:
					varLen = 0
				svTypes[record.info['SVTYPE']]={'count':1,'lengths':[varLen]}
			else:
				svTypes[record.info['SVTYPE']]={'count':1,'lengths':[record.info['SVLEN']]}
		else:
			svTypes[record.info['SVTYPE']]['count'] += 1
			if 'SVLEN' in record.info:
				svTypes[record.info['SVTYPE']]['lengths'].append(record.info['SVLEN'])
			else:
				ref_len = len(record.ref)
				alt_len = len(record.alts)
				if alt_len - ref_len == 0:
					varLen = ref_len
				else:
					varLen = 0
				svTypes[record.info['SVTYPE']]['lengths'].append(varLen)

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


	vcf_prefix = in_vcf.split(".")[0]
	print(f'finished analyzing VCF: {num_records} variants in the {vcf_prefix}' )
	for key, val in svTypes.items():
		print(key,val['count'])

	# write out the variant counts to a tsv file, including the header column names, but not the index
	sample_variant_count_df = pd.DataFrame( list(variant_counts.items()), columns=['Sample', 'VariantCount'])
	sample_variant_count_df = sample_variant_count_df.sort_values('VariantCount').reset_index(drop=True)
	
	sample_variant_count_df.to_csv(vcf_prefix+"_sample_variant_counts.tsv", header=True, index=False, sep="\t")


	# return the variant count dataframe to main
	return sample_variant_count_df, svTypes

def violin_swarm(x,y,data,ax,swarm_pt_size = 3):
	""" make a violin plot with a swarm of datapoints on top """

	sns.violinplot(x=x, y=y, data=data, cut=0.25, inner="quartile", alpha = 0.05, ax=ax, edgecolor='black')
	# sns.swarmplot(x=x, y=y, data=data, s=swarm_pt_size, alpha=1, ax=ax, color='black')
	sns.stripplot(
	    data=data,
	    x=x,
	    y=y,
	    jitter=True,   # random jitter
	    size=2, ax=ax, color='black'
	)

def violin_swarm_cohort(x,y,data,ax,swarm_pt_size = 3):
	""" make a violin plot with a swarm of datapoints on top """

	cohort_colors = {"PPMI" : '#2ca02c',
					 "RUSH" : '#d62728',
					 "HBCC" : '#ff7f0e',
					 "NABEC" : '#1f77b4'
					}
	sns.violinplot(x=x, y=y, data=data, cut=0.25, inner="quartile", alpha = 0.05, ax=ax, edgecolor='black')
	# sns.swarmplot(x=x, y=y, data=data, s=swarm_pt_size, alpha=1, ax=ax, color='black')
	sns.stripplot(
		data=data,
		x=x,
		y=y,
		hue=x,
		jitter=True,
		size=2,
		ax=ax,
		palette=cohort_colors,
		legend=False
	)


def plot_violin_perSample(vcf_data, vcf_prefix, plot_title):
	""" set up the figure to plot a violin """
	print('Plotting sample count violin')
	fig, axs = plt.subplots(figsize=(8,8))


	violin_swarm(['samples']*vcf_data.shape[0], 'VariantCount', vcf_data, axs)

	plt.xticks(fontsize=ticksize)
	plt.yticks(fontsize=ticksize)
	# axs.set_xlabel('sample',fontsize=labelsize)
	axs.set_ylabel('variantCount',fontsize=labelsize)
	fig.suptitle(f"SV Count Per Sample {plot_title}", fontsize=labelsize)

	plt.tight_layout()
	plt.savefig(vcf_prefix+"_sample_variant_counts.png", dpi=300)
	plt.close(fig)

def plot_violin_perSample_perCohort(vcf_data, vcf_prefix, plot_title):
	""" set up the figure to plot a violin """
	print('Plotting sample count violin')
	fig, axs = plt.subplots(figsize=(6,8))

	vcf_data['cohort'] = vcf_data['Sample'].str.split('_').str[0]

	violin_swarm_cohort('cohort', 'VariantCount', vcf_data, axs)

	plt.xticks(fontsize=12, rotation=45)
	plt.yticks(fontsize=14)
	# axs.set_xlabel('sample',fontsize=labelsize)
	axs.set_ylabel('variantCount',fontsize=labelsize)
	fig.suptitle(f"SV Count Per Sample By Cohort {plot_title}", fontsize=labelsize)

	plt.tight_layout()
	plt.savefig(vcf_prefix+"_cohort_sample_variant_counts.png", dpi=300)
	plt.close(fig)


def plot_violin_variantType(svTypes, vcf_prefix, plot_title):
	""" Plot variant type violin of variant lengths """
	print('Plotting variant type violin')

	# convert dictionary to long-form DF
	lfdata = []
	for svtype, vals in svTypes.items():
		if 'lengths' in vals.keys():
			if len(vals['lengths'])>1:
				for length in vals['lengths']:
					lfdata.append({'SVTYPE':svtype, 'SVLEN':length})
			else:
				lfdata.append({'SVTYPE':svtype, 'SVLEN':vals['lengths']})
	
	lfdf = pd.DataFrame(lfdata)
	lfdf = lfdf.sort_values('SVLEN').reset_index(drop=True)

	lfdf_10kmax = lfdf.loc[(lfdf['SVLEN']>=-10000) & (lfdf['SVLEN']<=10000)]
	lfdf_large = lfdf.loc[(lfdf['SVLEN']<-10000) | (lfdf['SVLEN']>10000)]

	lfdf.to_csv(vcf_prefix+"_variantType_counts.tsv", header=True, index=False, sep="\t")


	fig, axs = plt.subplots(1,2, figsize=(18,16))

	violin_swarm('SVTYPE', 'SVLEN', lfdf_10kmax, axs[0])
	violin_swarm('SVTYPE', 'SVLEN', lfdf_large, axs[1])

	for ax, subtitle in zip(axs, ["|SVLEN| \u2264 10kb", "|SVLEN| > 10kb"]):
		ax.set_title(subtitle, fontsize=labelsize)
		ax.set_xlabel("SVTYPE", fontsize=labelsize)
		ax.set_ylabel("SVLEN", fontsize=labelsize)
		ax.tick_params(axis='x', labelsize=ticksize)
		ax.tick_params(axis='y', labelsize=ticksize)

	fig.suptitle(f"SV Length Distribution per SV Type {plot_title}", fontsize=labelsize)

	plt.savefig(vcf_prefix+"_variant_counts_lengths.png", dpi=300)
	plt.close(fig)

# def plot_histogram_varLength( svLengths, vcf_prefix):
# 	""" plot histogram and then a scatter of variant length """




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

	parser.add_argument(
		"-m","--main_title",
		type=str,
		help="Title prefix for plots."
	)

	# add arugment for making a plot of sample variant counts 
	parser.add_argument(
		'--plot_violin_perSample', 
		action='store_true', 
		help='make a violin plot of the number of variants per sample'
	)

	# add arugment for making a plot of sample variant type counts 
	parser.add_argument(
		'--plot_violin_variantType', 
		action='store_true', 
		help='make a violin plot of the number of variants types for a single sample'
	)

	parser.add_argument(
		'--writeOutvariantTypes', 
		action='store_true', 
		help='write out variants types for a single sample'
	)



	if len(sys.argv) == 0:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	args = parser.parse_args()

	# Process the VCF file
	sample_variant_count_df, svTypes = vcfEntriesPerSample(args.in_vcf_file)

	#vcf prefix
	vcf_prefix = args.in_vcf_file.split(".")[0]
	plot_title = args.main_title

	if args.plot_violin_perSample:
		plot_violin_perSample(sample_variant_count_df, vcf_prefix, plot_title)

		plot_violin_perSample_perCohort(sample_variant_count_df, vcf_prefix, plot_title)

	if args.plot_violin_variantType:

		plot_violin_variantType(svTypes,vcf_prefix, plot_title)	

	if args.writeOutvariantTypes:
		svTypesDf = pd.DataFrame( [ (k,v['count']) for k,v in svTypes.items() ], columns=['type','count'] )
		totalVars = svTypesDf['count'].sum()
		svTypesDf.loc[len(svTypes)] = ['Total',totalVars]
		svTypesDf.to_csv(vcf_prefix+"_variantType_df.tsv", header=True, index=False, sep="\t")

	print('\nDone!\n')


