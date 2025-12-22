import argparse
import pandas as pd
import numpy as np
import pysam
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import datetime


# convert this to manhattan plot for QTL hits

labelsize=10
ticksize=9

sns.set(
    style="white",  # sets the background to white
    rc={
        "xtick.labelsize": labelsize,
        "ytick.labelsize": labelsize,
        "axes.labelsize": labelsize,
        "axes.titlesize": ticksize,
        "figure.facecolor": "white"  # ensures figure background is white
    }
)



def plot_barchart(df_long):
    
    

    # split up types
    panel1_types = ["Hifiasm_INS", "Hifiasm_DEL", "Shasta_INS", "Shasta_DEL"]
    df_p1 = df_long[df_long['svtype'].isin(panel1_types)]
    df_p2 = df_long[~df_long['svtype'].isin(panel1_types)]

    # plt.figure(figsize=(12,6))
    fig, axs = plt.subplots(1, 2, figsize=(18, 8)) 
    sns.barplot(data=df_p1, x='svtype', y='count', hue='sample',width=0.9, dodge=True, palette='Set2', ax = axs[0])

    axs[0].set_title("INS DEL Counts per Sample")
    axs[0].set_ylabel("Count")
    axs[0].set_xlabel("SV Type")
    axs[0].tick_params(axis='x', rotation=45)
    axs[0].get_legend().remove()

    sns.barplot(data=df_p2, x='svtype', y='count', hue='sample', palette='Set2', ax = axs[1])

    axs[1].set_title("SV Counts per Sample")
    axs[1].set_ylabel("Count")
    axs[1].set_xlabel("SV Type")
    axs[1].tick_params(axis='x', rotation=45)
    axs[1].get_legend().remove()

    # axs[0].set_xticks(range(len(df_p1['svtype'].unique())))
    axs[0].set_xticks([0,1,2,3])


    plt.tight_layout()
    plt.savefig(f"Shasta_Hifiasm_sampleSVs.png",dpi=300, facecolor='white', transparent=False)


def parse_input(file):

    asm = ''
    if file.startswith('shasta'):
        asm='Shasta_'
    else:
        asm='Hifiasm_'

    samples = {}
    current_sample = None
    current_counts = {}

    with open(file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # identify new data block
            if line.startswith("structuralVariants"):
                # save previous block
                if current_sample is not None:
                    samples[current_sample] = current_counts
                
                # new sample name
                current_sample = "_".join(line.split("/", 1)[1].split("_")[0:3])
                current_counts = {}
                continue
            
            # skip header
            if line.startswith("type"):
                continue
            
            # parse counts
            if "\t" in line:
                svtype, count = line.split("\t")
                svtype = svtype.replace("Total:", "Total")
                svtype = asm+svtype
                current_counts[svtype] = int(count)

    # save last block
    if current_sample is not None:
        samples[current_sample] = current_counts

    df = pd.DataFrame(samples).fillna(0)

    return df

if __name__ == "__main__":

    # Make an argument parser
    parser = argparse.ArgumentParser(description="Plot a manhttan of SVs and SNVs from a QTL.")

    # add argument for the input vcf file
    parser.add_argument(
        "-s","--shasta_count_tsv",
        type=str,
        required=True,
        help="Path to the input file to be plotted."
    )

    parser.add_argument(
        "-f","--hifiasm_count_tsv",
        type=str,
        required=True,
        help="Path to the input file to be plotted."
    )


    parser.add_argument(
        "-o","--output_directory",
        type=str,
        required=True,
        help="Path to the input qtl tsv to be analyzed."
    )

    # parser.add_argument(
    #     "-f","--pval_field",
    #     type=str,
    #     default="bh_fdr",
    #     help="Path to the input qtl tsv to be analyzed."
    # )

    # parser.add_argument(
    #     '--vertical', 
    #     action='store_true', 
    #     help='Turn manhttan vertical'
    # )


    args = parser.parse_args()

    if not os.path.isdir(args.output_directory):
        os.makedirs(args.output_directory, exist_ok=True)

    shasta = parse_input(args.shasta_count_tsv)
    # convert to long form for sns
    shasta_long = shasta.reset_index().melt(id_vars='index',
                                var_name='sample',
                                value_name='count')
    shasta_long.rename(columns={'index': 'svtype'}, inplace=True)
    shasta_long['assembler'] = 'Shasta'

    hifi = parse_input(args.hifiasm_count_tsv)
    hifi_long = hifi.reset_index().melt(id_vars='index',
                                var_name='sample',
                                value_name='count')
    hifi_long.rename(columns={'index': 'svtype'}, inplace=True)
    hifi_long['assembler'] = 'Hifiasm'

    df = pd.concat([shasta_long, hifi_long], ignore_index=True)


    plot_barchart(df)


    # qtl_df = pd.read_csv(args.in_qtl_tsv, sep="\t")

    # 'variant_id'
    # chr1_806320_806320_INS_102
    # chr1:986336:C:A
    # qtl_df[['chrom','start','end_ref','type_alt','size']] = qtl_df['variant_id'].str.split(r'[_:]', expand=True)  
    # qtl_df['log10p'] = -np.log10(qtl_df[args.pval_field])  

    # print(qtl_df.shape )
    # prefix = args.in_qtl_tsv.split(".")[0]
    # # manhattan_plot(qtl_df, args.cohort, args.output_directory, prefix, 'log10p', args.vertical )

    # plot_volcano(qtl_df, prefix, args.output_directory, 'log10p', alpha=0.3)