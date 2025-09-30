import io
import os
import sys
import subprocess
import argparse
import pandas as pd

"""
Script to move and rename data to Terra Data staging workspace

inputs: a TSV with samples as rows, first column is sample ID
		data column of files that are to be copied

example: python3 download_gs_data.py -i NABEC_cohort_harmonizedVCF_methyl_012025.tsv -c out_HarmPhase_methylationBed2_GRCh38
After running script: chmod u+x output.sh
./output.sh
files copied to output directory named for column
		
Author: Melissa Meredith
12/2024
"""


def download_data(tsv_path, download_column, use_sample_in_filename):
	print('tsv_path', tsv_path)
	print('downlad data from column:', download_column)

	if (not os.path.isfile(tsv_path)):
		print(f"Error: File '{tsv_path}' not found.")
		return


	try:

		print(f'making commands for: {download_column}')
		bash_script = open(f"{download_column}_gs_commands.sh", "w")

		data_df = pd.read_csv(tsv_path, sep="\t")
		# print(data_df.head())

		samplecol = list(data_df.columns)[0]
		print('sample column name', samplecol)

		if download_column not in list(data_df.columns):
			print(download_column, 'not in tsv')
			return


		# make commands for each sample 
		for idx, row in data_df.iterrows():
			
			sample_id = row[samplecol]

			gslink = row[download_column]
			# if link is empty skip it
			if pd.isna(gslink):
				continue

			if use_sample_in_filename:
				gsfilename = gslink.split("/")[-1]
				filename = f"{sample_id}_{gsfilename}"
				# make gs command with sample id filename
				gs_command = f"gsutil cp {gslink} ./{download_column}/{filename}"
			else:
				gs_command = f"gsutil cp {gslink} ./{download_column}/" 
			

			bash_script.write(f'{gs_command}\n\n')


		bash_script.close()


	except Exception as e:
		print(f'error reading tsv file: {tsv_path}')
		print(e)




if __name__ == "__main__":

    # Make an argument parser
    parser = argparse.ArgumentParser(description="Process a TSV file to copy and bgzip files from GCS.")
    parser.add_argument(
        "-i","--in_tsv_file",
        type=str,
        required=True,
        help="Path to the input gs link TSV file with header. The first column should be the sample IDs, subsequent columns are data to be copied."
    )

    parser.add_argument(
        "-c","--column_to_download",
        type=str,
        required=True,
        help="Column in tsv to make download commands."
    )

    parser.add_argument(
        "-s","--use_sample_in_filename",
        action="store_true",
        help="Use the sample id column to add sample id to the beginning of the downloaded filename."
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # Parse arguments
    args = parser.parse_args()

    # Process the TSV file
    download_data(args.in_tsv_file, args.column_to_download, args.use_sample_in_filename)
