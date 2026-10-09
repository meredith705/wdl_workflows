#!/usr/bin/python 

# filter_input_methylation_tables.py
# script to filter input methylation tables (extended average methylation per region over regions and samples) by samples and chromosomes, amongst other features

import argparse
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
import time
import psutil
import os
import re

# subroutine to parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Filter input methylation tables (extended average methylation per region over regions and samples) by samples and chromosomes, amongst other features.")
    # input methylation table
    parser.add_argument("--input_methylation_table", 
        required=True, help="Path to input methylation table (methylation per sample per aggregated region).")
    # output methylation table
    parser.add_argument("--output_methylation_table", 
        required=False, help="Path to output methylation table after filtering (methylation per sample per aggregated region).")
    # add reference name to find sample columns
    parser.add_argument("--reference_name",
        required=False,default="GRCh38",help="Reference name used to separate sample from non-sample columns in table.")
    # included chromosomes
    parser.add_argument("--included_chromosomes", 
        required=False, help="Chromosomes to include in output (list in text file, one chromosome per line).")
    # excluded chromosomes
    parser.add_argument("--excluded_chromosomes", 
        required=False, help="Chromosomes to exclude in output (list in text file, one chromosome per line).")
    # included samples
    parser.add_argument("--included_samples", required=False, help="Samples to include in output (list in text file, one chromosome per line).")
    # excluded samples
    parser.add_argument("--excluded_samples", 
        required=False, help="Samples to exclude in output (list in text file, one chromosome per line).")
    # print only samples and chromosomes in input table
    parser.add_argument("--print_input_samples_and_chromosomes_only", 
        required=False, action=argparse.BooleanOptionalAction, default=False, help="Print samples and chromosomes in input table to files (input_methylation_table.samples, .chromosomes) only and exit.")
    # also print samples and chromosomes in output table
    parser.add_argument("--print_output_samples_and_chromosomes", 
        required=False, action=argparse.BooleanOptionalAction, default=False, help="Print samples and chromosomes in output table to files (output_methylation_table.samples, .chromosomes) as well.")
    # impute mean values per region
    parser.add_argument("--impute_with_region_means", 
        required=False, action=argparse.BooleanOptionalAction, default=False, help="Impute N/A methylation values as mean per region (default false).")
    # filter out regions with more than 5% of samples missing
    parser.add_argument("--missing_methylation_rate_filter", 
        required=False, type=float, default=0.05, help="Remove region if more than this proportion of samples have no methylation calls; default 0.05.")
    # remove reference_modFraction from sample names
    parser.add_argument("--clean_sample_names", 
        required=False, action=argparse.BooleanOptionalAction, default=False, help="Output table with sample name only for sample columns (no _[reference]_modFraction).")
    # output in tensorQTL phenotype format
    parser.add_argument("--tensorQTL_format",
        required=False,action=argparse.BooleanOptionalAction,default=False,help="Output table in tensorQTL phenotype BED format (chromosome/start/end/phenotype ID/samples).")
    # return parsed arguments
    return parser.parse_args()
    
# main script subroutine
def main():
    # Parse input command line arguments
    args = parse_args()

    # load chromosomes to include to list from file
    if (args.included_chromosomes is not None):
        with open(args.included_chromosomes, 'r') as file:
            # Use a list comprehension to strip the newline character from each line
            include_chromosome_list = [line.strip() for line in file]

    # load chromosomes to exclude to list from file
    if (args.excluded_chromosomes is not None):
        with open(args.excluded_chromosomes, 'r') as file:
            # Use a list comprehension to strip the newline character from each line
            exclude_chromosome_list = [line.strip() for line in file]

    # load samples to include to list from file
    if (args.included_samples is not None):
        with open(args.included_samples, 'r') as file:
            # Use a list comprehension to strip the newline character from each line
            include_sample_list = [line.strip() for line in file]

    # load samples to exclude to list from file
    if (args.excluded_samples is not None):
        with open(args.excluded_samples, 'r') as file:
            # Use a list comprehension to strip the newline character from each line
            exclude_sample_list = [line.strip() for line in file]

    # load table
    input_meth_table_df=pd.read_csv(args.input_methylation_table,sep="\t")
    # get list of chromosomes included in table
    input_meth_chromosome_list=input_meth_table_df['#chrom'].unique()
    print('input table:\n', input_meth_table_df.head(30))

    # get list of samples included in table - note reference is defined to be GRCh38
    # find all columns with modFraction, then remove _GRCh38 and everything after
    input_meth_samples_list=[n.split('_'+args.reference_name,1)[0] for n in input_meth_table_df.columns[input_meth_table_df.columns.str.endswith("modFraction")]]

    # print samples and chromosomes included in table
    input_meth_sample_list_filepath=args.input_methylation_table + ".samples"
    input_meth_chromosome_list_filepath=args.input_methylation_table + ".chromosomes"

    # print sample list to file
    with open(input_meth_sample_list_filepath, 'w') as file:
        for item in input_meth_samples_list:
            file.write(str(item) + '\n')
    # print chromosome list to file
    with open(input_meth_chromosome_list_filepath, 'w') as file:
        for item in input_meth_chromosome_list:
            file.write(str(item) + '\n')
    # exit now if --print_input_samples_and_chromosomes_only set
    if (args.print_input_samples_and_chromosomes_only is True):
        # samples based on output sample list
        print("Total samples in input:",len(input_meth_samples_list))
        # regions based on row count
        print("Total regions in input:",input_meth_table_df.shape[0])
        # chromosomes based on chromosome list
        print("Total chromosomes in input:",len(input_meth_chromosome_list))
        print("Printed input table sample and chromosome lists to file - now exiting.")
        exit()

    # initialize output table
    output_meth_table_df=input_meth_table_df.copy()
    # filter on chromosomes
    if (args.included_chromosomes is not None):
        output_meth_table_df=output_meth_table_df[output_meth_table_df['#chrom'] in included_chromosome_list]
    if (args.excluded_chromosomes is not None):
        output_meth_table_df=output_meth_table_df[output_meth_table_df['#chrom'] not in excluded_chromosome_list]

    # filter on samples
    # first samples to include
    if (args.included_samples is not None):
        # separate table into samples and nonsamples column tables
        # nonsample columns do not have reference name in column name
        output_meth_table_nonsamples_df=output_meth_table_df[[item for item in output_meth_table_df.columns if args.reference_name not in item]]

        # sample columns have reference name in column name
        output_meth_table_samples_df=output_meth_table_df[[item for item in output_meth_table_df.columns if args.reference_name in item]]
        output_meth_table_sample_names=[n.split('_'+args.reference_name,1)[0] for n in output_meth_table_samples_df.columns]
        
        # get included sample indices
        included_sample_indices=[output_meth_table_sample_names.index(i) for i in include_sample_list]
        # only include these samples in samples df
        output_meth_table_samples_df=output_meth_table_samples_df.iloc[:,included_sample_indices]
        # concatenate samples df with nonsamples df
        output_meth_table_df=pd.concat([output_meth_table_nonsamples_df,output_meth_table_samples_df],axis=1)

    # then samples to exclude
    if (args.excluded_samples is not None):
        # separate table into samples and nonsamples column tables
        # nonsample columns do not have reference name in column name
        output_meth_table_nonsamples_df=output_meth_table_df[[item for item in output_meth_table_df.columns if args.reference_name not in item]]
        # sample columns have reference name in column name
        output_meth_table_samples_df=output_meth_table_df[[item for item in output_meth_table_df.columns if args.reference_name in item]]
        output_meth_table_sample_names=[n.split('_'+args.reference_name,1)[0] for n in output_meth_table_samples_df.columns]
        # get included sample indices
        excluded_sample_indices=[output_meth_table_sample_names.index(i) for i in exclude_sample_list]
        # only include these samples in samples df
        output_meth_table_samples_df=output_meth_table_samples_df.drop(columns=output_meth_table_samples_df.columns[excluded_sample_indices])
        # concatenate samples df with nonsamples df
        output_meth_table_df=pd.concat([output_meth_table_nonsamples_df,output_meth_table_samples_df],axis=1)

    # Replace `.` with `NaN`
    output_meth_table_df = output_meth_table_df.replace('.', np.nan)
    # filter regions on missing methylation rate after sample/chromosome filters
    if (args.missing_methylation_rate_filter is not None):
        # calculate missing methylation per row
        output_meth_table_samples_df=output_meth_table_df[[item for item in output_meth_table_df.columns if args.reference_name in item]]
        output_meth_table_samples_df=output_meth_table_samples_df.replace('.', np.nan).astype(float)
        output_meth_table_missing_methylation_info_rates=output_meth_table_samples_df.isna().sum(axis=1)/output_meth_table_samples_df.shape[1]
        print('missing ness rates:\n',output_meth_table_missing_methylation_info_rates)
        output_meth_table_df = output_meth_table_df[output_meth_table_missing_methylation_info_rates < args.missing_methylation_rate_filter]

    # impute missing values with means per region (after filtering on missing methylation rate!!)
    if (args.impute_with_region_means is True):
        # impute missing values per row with region means
        output_meth_table_samples_df=output_meth_table_df[[item for item in output_meth_table_df.columns if args.reference_name in item]]
        output_meth_table_samples_df=output_meth_table_samples_df.replace('.', np.nan).astype(float)
        output_meth_table_df_region_means = output_meth_table_samples_df.mean(axis=1)
        print('regional means:\n', output_meth_table_df_region_means.head(30))
        output_meth_table_df=output_meth_table_df.replace('.', np.nan)
        # output_meth_table_df.fillna(output_meth_table_df_region_means,axis=0)
        output_meth_table_df = output_meth_table_df.fillna(output_meth_table_df_region_means,axis=0)
        print('imputed table:\n', output_meth_table_df.head(15))

    # write output table samples and chromosomes to file if specified
    if (args.print_output_samples_and_chromosomes is True):
        # get list of chromosomes included in table
        output_meth_chromosome_list=output_meth_table_df['#chrom'].unique()
        # get list of samples included in table - note reference is defined to be GRCh38
        # find all columns with modFraction, then remove _GRCh38 and everything after
        output_meth_samples_list=[n.split('_'+args.reference_name,1)[0] for n in output_meth_table_df.columns[output_meth_table_df.columns.str.endswith("modFraction")]]
        # print samples and chromosomes included in table
        output_meth_sample_list_filepath=args.output_methylation_table + ".samples"
        output_meth_chromosome_list_filepath=args.output_methylation_table + ".chromosomes"
        # print sample list to file
        with open(output_meth_sample_list_filepath, 'w') as file:
            for item in output_meth_samples_list:
                file.write(str(item) + '\n')
        # print chromosome list to file
        with open(output_meth_chromosome_list_filepath, 'w') as file:
            for item in output_meth_chromosome_list:
                file.write(str(item) + '\n')
        # samples based on output sample list
        print("Total samples kept:",len(output_meth_samples_list))
        # regions based on row count - no wait, output chromosome list
        print("Total regions kept:",output_meth_table_df.shape[0])
        # chromosomes based on chromosome list
        print("Total chromosomes kept:",len(output_meth_chromosome_list))
    # reformat output BED to tensorQTL format if specified
    if (args.tensorQTL_format is True):
        print("Converting output table to tensorQTL standard BED format.")
        # first four columns are chromosome, start, end, and phenotype id
        # make phenotype id from chromosome, start, and end (check for uniqueness)
        # next are samples
        # change sample names to correct ones
        # print(output_meth_table_df)
        output_meth_table_tensorQTL_df=output_meth_table_df[[item for item in output_meth_table_df.columns if args.reference_name in item]]
        output_meth_table_sample_names=[n.split('_'+args.reference_name,1)[0] for n in output_meth_table_tensorQTL_df.columns]
        output_meth_table_tensorQTL_df.columns=output_meth_table_sample_names
        # print(output_meth_table_tensorQTL_df)
        #output_meth_table_tensorQTL_df['#chr']=output_meth_table_df['#chrom']
        output_meth_table_tensorQTL_df_bedcols=pd.DataFrame(columns=['#chr','start','end'])
        output_meth_table_tensorQTL_df_bedcols['#chr']=output_meth_table_df.iloc[:,0]
        # output_meth_table_tensorQTL_df['start']=output_meth_table_df['start']
        output_meth_table_tensorQTL_df_bedcols['start']=output_meth_table_df.iloc[:,1].astype(int)
        # output_meth_table_tensorQTL_df['end']=output_meth_table_df['end']
        output_meth_table_tensorQTL_df_bedcols['end']=output_meth_table_df.iloc[:,2].astype(int)
        # default phenotype id pattern is chrom_start_end - check for uniqueness after making table
        output_meth_table_tensorQTL_df_bedcols['phenotype_id']=output_meth_table_df.iloc[:,0].astype(str)+"_"+output_meth_table_df.iloc[:,1].astype(str)+"_"+output_meth_table_df.iloc[:,2].astype(str)
        # output_meth_table_tensorQTL_df['phenotype_id']=output_meth_table_df['#chrom']+"_"+output_meth_table_df['start']+"_"+output_meth_table_df['end']
        # set output_meth_table_df to output_meth_table_tensorQTL_df
        output_meth_table_df=pd.concat([output_meth_table_tensorQTL_df_bedcols,output_meth_table_tensorQTL_df],axis=1)
        # print(output_meth_table_df)
    # clean up sample names to sample name only, not with reference and modfraction suffix, if specified
    if (args.clean_sample_names is True):
        # nonsample columns do not have reference name in column name
        output_meth_table_nonsamples_df=output_meth_table_df[[item for item in output_meth_table_df.columns if args.reference_name not in item]]
        # sample columns have reference name in column name
        output_meth_table_samples_df=output_meth_table_df[[item for item in output_meth_table_df.columns if args.reference_name in item]]
        output_meth_table_sample_names=[n.split('_'+args.reference_name,1)[0] for n in output_meth_table_samples_df.columns]
        output_meth_table_samples_df.columns=output_meth_table_sample_names
        # concatenate samples df with nonsamples df
        output_meth_table_df=pd.concat([output_meth_table_nonsamples_df,output_meth_table_samples_df],axis=1)
    # print to stdout total number of samples and regions included in output
    print("Filtering complete, writing final output methylation table to file.")
    # write table to file
    output_meth_table_df.to_csv(args.output_methylation_table,sep="\t",index=False)

if __name__ == "__main__":
    main()