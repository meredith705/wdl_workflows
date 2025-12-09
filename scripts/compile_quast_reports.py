#!/usr/bin/env python3

import os
import argparse
import pandas as pd

def parse_quast_report(report_path):
    """
    Parse a QUAST report.txt file into a dictionary of metrics.
    """
    metrics = {}
    with open(report_path) as f:
        for line in f:
            # skip first 3 lines in report.txt
            if line.strip() == "" or line.startswith("A"):
                continue
            # Each line looks like: "Genome fraction (%)   97.5"
            # split on whitespace
            parts = line.strip().split()

            if len(parts)>2:
                asm_names = ["hapdup_dual_1", "hapdup_dual_2"]
            else:
                asm_names = ["hapdup_dual_1"]

            # quast can have 2 assembly columns now 
            key = " ".join(parts[:-2]) if len(parts) > 2 else parts[0]
            #key = " ".join(parts[:-1])
            
            # store values for 2 assemblies
            if len(asm_names) ==1:
                value = parts[-1] 
            else: 
                value1 = parts[-2]
                value2 = parts[-1]

            # # unaligned contigs is 1 + 577 part
            if ("# unaligned contigs" in line) or ( "# genomic features" in line ):
                print(line)
                # store key as the first space separated parts of the line 
                key = parts[0] + " " + parts[1] + " " + parts[2]

                # for quast reports with just one assembly store one value set
                if len(asm_names) ==1:
                    value = " ".join(parts[-4:])
                else:
                    # for 2 assemblies store asm1 and asm2 value sets
                    value1 = " ".join(parts[3:-4])
                    value2 = " ".join(parts[-4:])
           
            if len(asm_names) ==1:
                metrics[key] = value
            else: 
                metrics[key+"_"+asm_names[0]] = value1
                metrics[key+"_"+asm_names[1]] = value2

    return metrics

def main(report_dir, output_file):
    """
    Collect metrics from multiple QUAST reports and compile into a table.
    """
    all_reports = {}
    for root, dirs, files in os.walk(report_dir):
        for fname in files:
            file_parts = fname.split("_")
            sample_name = "_".join(file_parts[:-1])
            file_name = file_parts[-1]
            # check the filename and store metrics in dictionary by sample name
            if file_name == "report.txt":
                report_path = os.path.join(root, fname)
                print('report_path', report_path, fname, sample_name)
                metrics = parse_quast_report(report_path)
                all_reports[sample_name] = metrics

    # Convert to dataframe
    df = pd.DataFrame.from_dict(all_reports, orient="index")
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
