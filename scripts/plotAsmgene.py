#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

"""
Plots asmgene data for SHasta Hifiasm compariosn. 

First run: 
python3 ~/Desktop/paten/wdl_git_workflows/wdl_workflows/scripts/download_gs_data.py -i ../../ROSMAP_cohort.tsv -c AsmgenePerGeneStats --use_sample_in_filename

chmod u+x AsmgeneStats_gs_commands.sh
./AsmgeneStats_gs_commands.sh

python3 ~/Desktop/paten/wdl_git_workflows/wdl_workflows/scripts/compile_asmgene_report.py -i AsmgeneStats -o AsmgeneStats_RUSH_CHM13_0.97.compiledASMGene.tsv

python3 ~/Desktop/paten/wdl_git_workflows/wdl_workflows/scripts/plotAsmgene.py 
    -i ShastaHD_AsmgeneStats_RUSH_CHM13_0.97.compiledASMGene.tsv Hifiasm_AsmgeneStatsHifiasm_RUSH_CHM13_0.97.compiledASMGene.tsv -o asmgene_shastaHD_hifiasm_MMC_PSC.png

"""


def main(input_tsvs, output_png, plot_sampleNames, customGroups):
    plt.figure(figsize=(8, 6))
    # fig, axes = plt.subplots(1, 2, figsize=(10, 8))
    plt.style.use('ggplot') 

    assembler_colors = {
        "ShastaHD": "steelblue",
        "Hifiasm": "darkorange"
    }

    if customGroups:
        # Get the colors from the current style's property cycle
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        if len(customGroups) > len(colors):
            print("not enough colors:", len(colors))
            return 

        # for i, group in enumerate(customGroups):
        assembler_colors = {
                customGroups[i]: colors[i] for i in range(len(customGroups))
            }


    for i, tsv in enumerate(input_tsvs):
        # Load compiled QUAST report
        df = pd.read_csv(tsv, sep="\t")
        df = df.rename(columns={'Unnamed: 0': 'Sample'})

        if customGroups:
            assembler = customGroups[i]
        else:
            assembler = tsv.split("_")[0]
        
        print('assembler plotting:', assembler)

        df['Assembler'] = assembler

        # Clean up column names (QUAST often includes spaces)
        # df.columns = df.columns.str.strip()

        # # Ensure numeric types
        # df["# misassemblies_hapdup_dual_1"] = pd.to_numeric(df["# misassemblies_hapdup_dual_1"], errors="coerce")
        # df["# mismatches per 100 kbp_hapdup_dual_1"] = pd.to_numeric(df["# mismatches per 100 kbp_hapdup_dual_1"], errors="coerce")

        # Scatterplot: mismatches (x) vs misassemblies (y), hue = assembler, style = assembler
        
        sns.scatterplot(
            data=df,
            x="pctMMC",
            y="pctSingleCopy",
            hue="Assembler",
            style="Assembler",
            s=80, alpha=0.95,
            palette=assembler_colors
        )

        if plot_sampleNames:
            # Annotate sample names
            for _, row in df.iterrows():
                plt.text(
                    row["pctMMC"] + 0.001,
                    row["pctSingleCopy"] + 0.001,
                    row["Sample"],
                    fontsize=8
                )

    plt.xlabel("Percent Missing Multicopy Genes")
    plt.ylabel("Percent of Single Copy Genes Assembled")
    plt.title("Asmgene Transcripts Mapped to CHM13")
    plt.legend(title="Assembler")
    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    plt.close()
    print(f"Plot saved to {output_png}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scatterplot of QUAST misassemblies vs mismatches per 100 kbp")
   
    parser.add_argument(
        "-i", "--input_tsvs",
        nargs="+",
        required=True,
        help="List of compiled Asmgene TSV files"
    )

    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Output plot image (PNG)")

    parser.add_argument(
        '--plot_sampleNames', 
        action='store_true', 
        help='Plot sample names by point'
    )

    parser.add_argument(
        '--groups',
        nargs='+',
        help='A list of assembly types'
    )

    args = parser.parse_args()
    main(args.input_tsvs, args.output, args.plot_sampleNames, args.groups) 

