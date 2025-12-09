#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

"""
Plots quast data for SHasta Hifiasm compariosn. 

First run: 
python3 ~/Desktop/paten/wdl_git_workflows/wdl_workflows/scripts/download_gs_data.py -i ROSMAP_cohort.tsv -c quast_asmDual1_2_Q10_R10_k24_CHM13 --use_sample_in_filename

./quast_asmDual1_2_Q10_R10_k24_CHM13_gs_commands.sh

python3 ~/Desktop/paten/wdl_git_workflows/wdl_workflows/scripts/compile_quast_reports.py -i quast_asmDual1_2_Q10_R10_k24_CHM13 -o quast_asmDual1_2_Q10_R10_k24_CHM13.compiledQUAST.tsv

"""


def main(input_tsvs, output_png):
    # plt.figure(figsize=(8, 6))
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    assembler_colors = {
        "shasta": "steelblue",
        "hifiasm": "darkorange"
    }

    for tsv in input_tsvs:
        # Load compiled QUAST report
        df = pd.read_csv(tsv, sep="\t")
        df = df.rename(columns={'Unnamed: 0': 'Sample'})

        assembler = tsv.split("_")[0]
        print('assembler plotting:', assembler)

        df['Assembler'] = assembler

        # Clean up column names (QUAST often includes spaces)
        df.columns = df.columns.str.strip()

        # Ensure numeric types
        df["# misassemblies_hapdup_dual_1"] = pd.to_numeric(df["# misassemblies_hapdup_dual_1"], errors="coerce")
        df["# mismatches per 100 kbp_hapdup_dual_1"] = pd.to_numeric(df["# mismatches per 100 kbp_hapdup_dual_1"], errors="coerce")

        # Scatterplot: mismatches (x) vs misassemblies (y), hue = assembler, style = assembler
        
        sns.scatterplot(
            data=df,
            x="# mismatches per 100 kbp_hapdup_dual_1",
            y="# misassemblies_hapdup_dual_1",
            hue="Assembler",
            style="Assembler",
            s=100,
            palette=assembler_colors
        )

        # Annotate sample names
        for _, row in df.iterrows():
            plt.text(
                row["# mismatches per 100 kbp_hapdup_dual_1"] + 0.1,
                row["# misassemblies_hapdup_dual_1"] + 0.1,
                row["Sample"],
                fontsize=8
            )

    plt.xlabel("Mismatches per 100 kbp")
    plt.ylabel("# Misassemblies")
    plt.title("QUAST comparison: Misassemblies vs Mismatches per 100 kbp")
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
        help="List of compiled QUAST TSV files"
    )


    parser.add_argument(
        "-o", "--output", 
        required=True, 
        help="Output plot image (PNG)")

    args = parser.parse_args()
    main(args.input_tsvs, args.output)
