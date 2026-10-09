#!/usr/bin/env python3

import os
import argparse
import pandas as pd


HUMAN_GENOME_GB = 3.1  # human genome ~3.1 Gbp

def main(report_dir, output_file):
    """
    Collect metrics from multiple wambam-summary.csv reports and compile into a table.
    """
    all_reports = {}
    for root, dirs, files in os.walk(report_dir):
        for fname in files:
            file_parts = fname.split("_")
            sample_name = "_".join(file_parts[:-1])
            file_name = file_parts[-1]
            # check the filename and store metrics in dictionary by sample name
            if file_name == "wambam-summary.csv":
                report_path = os.path.join(root, fname)
                print('report_path', report_path, fname, sample_name)

                metrics = pd.read_csv(report_path)
                metrics.index = [sample_name]
                metrics['coverage'] = metrics['total.gbp']/HUMAN_GENOME_GB
                # print(metrics.head())
                all_reports[sample_name] = metrics

    # Convert to dataframe
    df = pd.concat(all_reports.values()).sort_index()
    print(df.head())

    # Save to file
    df.to_csv(output_file, sep="\t")
    print(f"Compiled report written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compile wambam-summary.csv files")
    parser.add_argument(
        "-i", "--input_dir",
        required=True,
        help="Directory containing wambam-summary outputs (subdirectories with <sample>_wambam-summary.csv)"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file (TSV format)"
    )
    args = parser.parse_args()

    main(args.input_dir, args.output)
