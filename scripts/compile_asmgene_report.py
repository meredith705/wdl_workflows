#!/usr/bin/env python3

import os
import argparse
import pandas as pd

"""
First run: 
python3 ~/Desktop/paten/wdl_git_workflows/wdl_workflows/scripts/download_gs_data.py -i ../../ROSMAP_cohort.tsv -c AsmgenePerGeneStats --use_sample_in_filename

chmod u+x AsmgeneStats_gs_commands.sh
./AsmgeneStats_gs_commands.sh

python3 ~/Desktop/paten/wdl_git_workflows/wdl_workflows/scripts/compile_asmgene_report.py -i AsmgeneStats -o AsmgeneStats_RUSH_CHM13_0.97.compiledASMGene.tsv

"""

def parse_asmgene_gene_stats(report_path):
    """
    Parse a asmgene gene_stats.txt file into a dictionary of metrics.
    """
    metrics = {}
    with open(report_path) as f:
        for line in f:
            # skip first line, though it has headers
            if line.strip() == "":
                continue

            headers = ['H', 'Metric', 'genesToRef', 'hapdup_dual_1']
            if line.startswith("H"):
                headers = line.strip().split()

            # Each line looks like: "X  full_sgl    34173   31413"
            # split on whitespace
            parts = line.strip().split()


            # asmgene has 4 columns  
            # key = " ".join(parts[:-2]) if len(parts) > 2 else parts[0]
            key = parts[1]

            # ref value
            ref_value = parts[2]
            # asm value
            asm_value = parts[3]
                   
            # ref metrics
            metrics['ref_'+key] = ref_value
            metrics['asm_'+key] = asm_value


    return metrics

def main(report_dir, output_file):
    """
    Collect metrics from multiple QUAST reports and compile into a table.
    """
    all_reports = {}
    for root, dirs, files in os.walk(report_dir):
        for fname in files:
            file_parts = fname.split("_")
            sample_name = "_".join(file_parts[:3])
            file_name = file_parts[-2:]
            # check the filename and store metrics in dictionary by sample name
            if fname.endswith("gene_stats.txt"):
                report_path = os.path.join(root, fname)
                print('report_path', report_path, fname, sample_name)
                metrics = parse_asmgene_gene_stats(report_path)
                all_reports[sample_name] = metrics

    # Convert to dataframe
    df = pd.DataFrame.from_dict(all_reports, orient="index")
    df['asm_dup_cnt'] = pd.to_numeric(df['asm_dup_cnt'], errors='coerce')
    df['ref_dup_cnt'] = pd.to_numeric(df['ref_dup_cnt'], errors='coerce')
    df['asm_full_sgl'] = pd.to_numeric(df['asm_full_sgl'], errors='coerce')
    df['ref_full_sgl'] = pd.to_numeric(df['ref_full_sgl'], errors='coerce')

    # compute % missing multi-copy genes
    df['pctMMC'] = 1 - (df['asm_dup_cnt']/df['ref_dup_cnt'])
    df['pctSingleCopy'] = df['asm_full_sgl']/df['ref_full_sgl']
    print(df.head())

    # Save to file
    df.to_csv(output_file, sep="\t")
    print(f"Compiled report written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compile QUAST report.txt files")
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Directory containing QUAST outputs (subdirectories with report.txt)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file (TSV format)"
    )
    args = parser.parse_args()

    main(args.input_dir, args.output)
