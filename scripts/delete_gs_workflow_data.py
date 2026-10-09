import io
import os
import sys
import subprocess
from pathlib import PurePosixPath
import argparse
import pandas as pd

"""
Script to delete Terra Data 

inputs: a google workspace id, 
		a workflow submission id, 
		the name of a workflow whos output is to be deleted, 
		the file ending of files to be deleted, 


example: python3 delete_gs_workflow_data.py -w <gs://fc-secure..> -s <submission_id> -k <workflow_id> -e .vcf.gz


After running script: chmod u+x *_gs_delete_commands.sh
confirm the correct files are listed for deletion.
and run the commands:
./*_gs_delete_commands.sh


Author: Melissa Meredith
1/2026
"""


def delete_data(tsv_path, delete_column):
	print('tsv_path', tsv_path)
	print('delete data from column:', delete_column)

	if not os.path.isfile(tsv_path):
		print(f"Error: File '{tsv_path}' not found.")
		return


	try:

		print(f'making commands for: {delete_column}')
		bash_script = open(f"{delete_column}_gs_commands.sh", "w")

		data_df = pd.read_csv(tsv_path, sep="\t")
		# print(data_df.head())

		samplecol = list(data_df.columns)[0]
		print('sample column name', samplecol)

		if delete_column not in list(data_df.columns):
			print(delete_column, 'not in tsv')
			return


		# make commands for each sample 
		for idx, row in data_df.iterrows():
			
			sample_id = row[samplecol]

			gslink = row[delete_column]

			# make gs command 
			gs_command = f"gsutil rm {gslink} "

			bash_script.write(f'{gs_command}\n\n')


		bash_script.close()


	except Exception as e:
		print(f'error reading tsv file: {tsv_path}')
		print(e)


def delete_directory_data(gs_workspace, submission_id, workflow, fileEnding):
	''' delete data from intermediate directories that 
		remained after a workflow 
		ex: delete_directory_data(###, indexSingleBam, sorted.bam)
		
		Given a submission id, list all submission/exection subdirectories, 
		store gs link to all files with provided file ending ( or some other
		identifier )

		Output a script of delete commands for review by user, then 
		run that bash script to delete data. 
		
		gs://fc-secure-dfed8d18-05dd-451e-ba76-75e17c3fc7e9/submissions/7fb66ca0-2258-4e9e-bd96-f9ede2938925/
		cardEndToEndVcfMethyl/b8184811-6b37-4cb9-b2c5-ec880cc5c05c/call-indexSingleInputBam/
		RUSH_007_FTX.sorted_minimap2_bq10_filtered.sorted.bam
	
		RUSH workspace ID: gs://fc-secure-dfed8d18-05dd-451e-ba76-75e17c3fc7e9
		ex submission id: 7fb66ca0-2258-4e9e-bd96-f9ede2938925
	'''

	print(f'delete {fileEnding} files from {workflow} in submission_id {submission_id} of {gs_workspace}')

	# build the gs link to the submission with files to be deleted
	baselocation = f"{gs_workspace}/submissions/{submission_id}"

	# make change file ending for output file nomenclature
	stringEnding = fileEnding.lstrip(".").replace(".","_")
	bash_script_outfile = f"{workflow}_{stringEnding}_gs_delete_commands.sh"
	print(f'writing delete commands to:\n{bash_script_outfile}')

	# stream output to save on memory
	proc = subprocess.Popen(
					["gsutil", "ls", "-r", baselocation],
					stdout=subprocess.PIPE,
					stderr=subprocess.PIPE,
					text=True
				)


	# isolate links to files that are to be deleted and write to file
	with open(bash_script_outfile, "w") as bash_script:

		for line in proc.stdout:

			gspath = line.strip() #.split("/")

			if not gspath.endswith(fileEnding):
				continue

			pathParts = PurePosixPath(gspath).parts
			

			if workflow in pathParts:
				print('pathParts', pathParts)
				# make gs command 
				gs_command = f"gsutil rm {gspath} "

				bash_script.write(f'{gs_command}\n\n')


	proc.stdout.close()
	rc = proc.wait()
	if rc != 0:
		raise RuntimeError("gsutil ls -r failed:\n{stderr}")

	print(f'finished writing commands.')


if __name__ == "__main__":
	"""
	delete_gs_workflow_data.py -w <gs://fc-secure..> -s <submission_id> -k <workflow_id> -e .vcf.gz
	"""

	# Make an argument parser
	parser = argparse.ArgumentParser(description="Produce commands to delete files in a google workspace by submission, workflow id and file ending.")
	
	parser.add_argument(
		"-w","--workspace_id",
		type=str,
		required=True,
		help="google workspace id <gs://...>."
	)

	parser.add_argument(
		"-s","--submission_id",
		type=str,
		required=True,
		help="Submission id of the main workflow."
	)

	parser.add_argument(
		"-k","--workflow_id",
		type=str,
		required=True,
		help="Workflow id of the wdl workflow with files to be deleted."
	)

	parser.add_argument(
		"-e","--file_ending",
		type=str,
		required=True,
		help="File ending of the desired files to be deleted. (eg. .bam, .vcf.gz)"
	)

	if len(sys.argv) == 1:
		parser.print_help(sys.stderr)
		sys.exit(1)

	# Parse arguments
	args = parser.parse_args()

	# Process the workflow submission data
	delete_directory_data(args.workspace_id, args.submission_id, 
							args.workflow_id, args.file_ending)

