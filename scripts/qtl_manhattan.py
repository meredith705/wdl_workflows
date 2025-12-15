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

labelsize=16
ticksize=14

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

def manhattan_plot(df, cohort, directory_path, prefix, pval, vertical):   
    
    #pval = 'bh_fdr'
        

    # List of 'chrom' values to color
    chrom_colors = {
        'chr1': 'red',
        'chr2': 'green',
        'chr3': 'blue',
        'chr4': 'orange',
        'chr5': 'purple',
        'chr6': 'brown',
        'chr7': 'pink',
        'chr8': 'gray',
        'chr9': 'cyan',
        'chr10': 'yellow',
        'chr11': 'magenta',
        'chr12': 'olive',
        'chr13': 'lime',
        'chr14': 'teal',
        'chr15': 'blueviolet',
        'chr16': 'gold',
        'chr17': 'navy',
        'chr18': 'indigo',
        'chr19': 'darkred',
        'chr20': 'darkgreen',
        'chr21': 'darkblue',
        'chr22': 'darkorange',
        'chrX': 'darkcyan',
        'chrY': 'orangered',
        'chrM': 'black'
    }

    colors=['#4D4D4D', '#B3B3B3']
    
    # make numeric chrom names
    # df['numeric_chrom'] = df.index.get_level_values('chrom').str.strip("chr")
    df['numeric_chrom'] = df['chrom'].str.strip("chr")
    df['numeric_chrom'] = pd.to_numeric(df['numeric_chrom'], errors='coerce')
    df_sorted = df.sort_values(by='numeric_chrom', ascending=False)
    # df_sorted.to_csv(directory_path+"/"+cohort.upper()+"_"+directory_path+"_variance_std.csv", header=True, index=True, sep=',')


    # fig, (ax_combined, axs, ax_std_deviation) = plt.subplots(3, 1, figsize=(10, 12), sharex=False)
    if vertical:
        w = 4
        h = 12
    else:
        w=10
        h = 5

    fig, axs = plt.subplots(figsize=(w,h))


    # Scatter plot for qtl ( manhattan )
#     axs.scatter(np.arange(variance_per_column.shape[0]), variance_per_column, color='blue', alpha=0.7, s=2)
    xtix = []
    xtix_labels = []
    last_group_end = 0
    colori = 0
    
    for chrom, group in df_sorted.groupby('numeric_chrom'):
        # axs.scatter(group['start'], group[pval], color=chrom_colors[chrom], alpha=0.4, s=2)
        if vertical:
            axs.scatter( group[pval], group['start'], color=colors[colori%2], alpha=1, s=16, label='SNV')
            xtix.append( (last_group_end+axs.get_yticks()[-2])/2)
            xtix_labels.append( 'chr'+str(chrom) )
            last_group_end = axs.get_yticks()[-1]
        else:
            axs.scatter(group['start'], group[pval], color=colors[colori%2], alpha=1, s=2, label='SNV')
            xtix.append( (last_group_end+axs.get_xticks()[-2])/2)
            xtix_labels.append( 'chr'+str(chrom) )
            last_group_end = axs.get_xticks()[-1]
        
        colori+=1

    # filtered = df_sorted.loc[df_sorted['type_alt'].isin(['INS', 'DEL'])]
    filtered = df_sorted.loc[df_sorted['CAVIAR Top Variant Type (SV>=SNV probability)']=='SV']
    for chrom, group in filtered.groupby('numeric_chrom'):
        if vertical:
            axs.scatter(group[pval], group['start'], color='red', alpha=0.7, s=20, label='SV')
        else:    
            axs.scatter(group['start'], group[pval], color='red', alpha=0.7, s=20, label='SV')
        # xtix.append( (last_group_end+axs.get_xticks()[-2])/2)
        # xtix_labels.append( chrom )
        # last_group_end = axs.get_xticks()[-1]
        
    # median_p = df[pval].median()
    median_p = -np.log10(0.05)
    if vertical:
        # axs.axvline(median_p, color='darkblue', alpha=0.3, linestyle='--', label=f'0.05 BH-FDR p-value')
        axs.set_xlabel('-log10(p)', fontsize = labelsize)
        axs.set_yticks(xtix)
        axs.set_yticklabels(xtix_labels)
    else:
        axs.axhline(median_p, color='darkblue', alpha=0.3, linestyle='--', label=f'0.05 BH-FDR p-value')
        axs.set_ylabel('-log10(p)', fontsize = labelsize)
        axs.set_xticks(xtix)
        axs.set_xticklabels(xtix_labels,rotation=90)
    # axs.legend(prop ={'size': labelsize})
    legend_element = [Line2D([0], [0], marker='o', color='red', label='SV', 
                         markerfacecolor='red', markersize=8,linestyle='None'),
                      Line2D([0], [0], marker='o', color=colors[0], label='SNV', 
                         markerfacecolor=colors[0], markersize=8,linestyle='None'),
                      Line2D([0], [0], marker='o', color=colors[1], label='SNV', 
                         markerfacecolor=colors[1], markersize=8,linestyle='None'),] 

    plt.legend(handles=legend_element)

    plt.gca().invert_yaxis()

    # Set title for the combined plot
    axs.set_title(f'{" ".join(prefix.split("_")[1::])}', fontsize=18)

    # Show the plots
    plt.tight_layout()
    
    plt.savefig(directory_path+"/"+prefix+"_"+directory_path+"_qtl_manhattan.png",dpi=300)

def plot_volcano(df, prefix, directory_path, pval, alpha=0.3):
    """ Linear Regression Volcano Plot""" 

    print(f"Plotting volcano")

    fig, axs = plt.subplots(1, 1, figsize=(10, 10)) 

    df['size'] = df['size'].fillna(1)
    df['size'] = pd.to_numeric(df['size'], errors='coerce')
    df_sorted = df.sort_values(by='size', ascending=True) 

    sns.scatterplot(x='slope', y=pval, data=df_sorted, color='#B3B3B3', s=16, ax=axs)
    sns.scatterplot(x='slope', y=pval, data=df_sorted.loc[(df_sorted['size']>20) | (df_sorted['size']<-20)], hue='size', edgecolor='k', alpha=0.7, s=18, ax=axs, palette='viridis')
    # bh_num_sig = df.loc[(df['bh_corrected_P_Value_Age']<alpha)].shape[0]
    axs.axhline(y=-np.log10(alpha), color='red', label=f'FDR={alpha}')

    # axs.legend()
    axs.set_title(f'{" ".join(prefix.split("_")[1::])}', fontsize=18)
    axs.set_ylabel("B-H Corrected - log10(p)")
    axs.set_xlabel("Slope")
    axs.grid(True)

    # t = plt.xticks(rotation=90)

    plt.savefig(f"{directory_path}/{prefix}_volcano.png",dpi=300, facecolor='white', transparent=False)


if __name__ == "__main__":

    # Make an argument parser
    parser = argparse.ArgumentParser(description="Plot a manhttan of SVs and SNVs from a QTL.")

    # add argument for the input vcf file
    parser.add_argument(
        "-i","--in_qtl_tsv",
        type=str,
        required=True,
        help="Path to the input qtl tsv to be analyzed."
    )

    parser.add_argument(
        "-c","--cohort",
        type=str,
        required=True,
        help="cohort ID."
    )

    parser.add_argument(
        "-o","--output_directory",
        type=str,
        required=True,
        help="Path to the input qtl tsv to be analyzed."
    )

    parser.add_argument(
        "-f","--pval_field",
        type=str,
        default="bh_fdr",
        help="Path to the input qtl tsv to be analyzed."
    )

    parser.add_argument(
        '--vertical', 
        action='store_true', 
        help='Turn manhttan vertical'
    )


    args = parser.parse_args()

    if not os.path.isdir(args.output_directory):
        os.makedirs(args.output_directory, exist_ok=True)

    qtl_df = pd.read_csv(args.in_qtl_tsv, sep="\t")

    # 'variant_id'
    # chr1_806320_806320_INS_102
    # chr1:986336:C:A
    qtl_df[['chrom','start','end_ref','type_alt','size']] = qtl_df['variant_id'].str.split(r'[_:]', expand=True)  
    qtl_df['log10p'] = -np.log10(qtl_df[args.pval_field])  

    print(qtl_df.shape )
    prefix = args.in_qtl_tsv.split(".")[0]
    # manhattan_plot(qtl_df, args.cohort, args.output_directory, prefix, 'log10p', args.vertical )

    plot_volcano(qtl_df, prefix, args.output_directory, 'log10p', alpha=0.3)










