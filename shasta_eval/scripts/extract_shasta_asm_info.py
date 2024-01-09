#!/usr/bin/env python3
"""
Author: Melissa Meredith
10/11/2023
A Script to QC Shasta Assemblies by extracting assembly and read stats reported in the shasta.log.tar.gz
stores as part of the NAPU Shasta.wdl or the NAPU end2end.wdl.

Can be run on one shasta.log.tar.gz file or a list of file paths.

Example usage to find file paths and run the script on that file path

find . -maxdepth 2 -not -path '*/\.*' | sed 's/^\.\///g' | grep "shasta.log.tar.gz" >> shataLogFiles.txt

python3.9 extract_shasta_asm_info.py -f shataLogFiles.txt

"""

import argparse
import os
import tarfile
import matplotlib.pyplot as plt


# make lists to store values
samples = []
read_n50s = []
shasta_n50s = []
shasta_assembled_len = []
estimated_genome_coverage = []


def get_stasta_asm_info(infile):
    """ Extract the read and assembly info reported by Shasta """
    added_read_n50 = False
    added_gcov = False
    with open(infile, 'r') as file:
        for line in file.readlines():
            # isolate input file/sample name
            if "shasta --input" in line:
                tokens = line.strip().split()
                samp = tokens[2].split("/")[-1].split("_")[:3]
                samples.append("_".join(samp))

            # Store genome coverage as Gibase
            if "raw bases" in line:
                tokens = line.strip().split()
                if tokens[-1][:-1].isdigit():
                    if not added_gcov:
                        estimated_genome_coverage.append(int(tokens[-1][:-1])/3200000000)
                        added_gcov = True
                    else:
                        estimated_genome_coverage[-1]=int(tokens[-1][:-1])/3200000000
            # store read and assembly N50
            if "N50" in line:
                tokens = line.strip().split()
                if tokens[2] == 'read' and tokens[5].isdigit():
                    if not added_read_n50:
                        read_n50s.append(int(tokens[5]))
                        added_read_n50=True
                    else:
                        read_n50s[-1]=int(tokens[5])
                if tokens[2] == 'assembly' and tokens[5].isdigit():
                    shasta_n50s.append(int(tokens[5]))
            # store total assembly length
            if "Total length" in line:
                tokens = line.strip().split()
                if tokens[-1].isdigit():
                    # print('shasta assembled sequence', tokens[-1])
                    shasta_assembled_len.append(int(tokens[-1]))

def best_worst_labels():
    labels = [""]*len(shasta_n50s)
    worst_quarter = sorted(shasta_n50s)[:round(len(shasta_n50s)*.25)]
    best_quarter = sorted(shasta_n50s)[round(len(shasta_n50s) * .75):]
    if len(worst_quarter)==0:
        worst_quarter=[sorted(shasta_n50s)[0]]
    if len(best_quarter)==0:
        best_quarter=[sorted(shasta_n50s)[-1]]
    # print(round(len(shasta_n50s) * .75),best_quarter, worst_quarter, '\n', set(worst_quarter+best_quarter))
    # shasta_n50s[:worst_quarter])
    for i,n50 in enumerate(shasta_n50s):
        if n50 in set(worst_quarter+best_quarter):
            labels[i]=samples[i]
    return labels


def plotting(output_dir, add_lables=True, label_cutoff=.05):
    if add_lables:
        labels = best_worst_labels()
    else:
        labels = [""]*len(shasta_n50s)

    """ Plot a 2x2 scatterplot figure of Asm N50, Read N50, Coverage, Asm Length """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 10))#, sharex=True, sharey=True)
    # plot all the data we extracted from shasta logs
    ax1.plot(read_n50s, shasta_n50s, "bo")
    ax1.set(ylabel="Haploid Shasta N50",xlabel="Read N50", title="Read N50 by Haploid Shasta N50")
    ax1_txt =  [ax1.annotate(txt, (read_n50s[i], shasta_n50s[i]), xytext=(read_n50s[i]+0.2,
                                                shasta_n50s[i]+0.2)) for i, txt in enumerate(labels)]

    ax2.plot(estimated_genome_coverage, shasta_n50s, "bo")
    ax2.set(ylabel="Haploid Shasta N50", xlabel="estimated read genome coverage",
            title="Estimated Read Genome Cov by Haploid Shasta N50")
    [ax2.annotate(txt, (estimated_genome_coverage[i], shasta_n50s[i]), xytext=(estimated_genome_coverage[i]+0.2,
                                                shasta_n50s[i]+0.2)) for i, txt in enumerate(labels)]

    ax3.plot(read_n50s, estimated_genome_coverage, "bo")
    ax3.set(ylabel="estimated read genome coverage", xlabel="read N50",
            title="Read N50 by Haploid Shasta length")
    [ax3.annotate(txt, (read_n50s[i], estimated_genome_coverage[i]), xytext=(read_n50s[i]+0.2,
                                                estimated_genome_coverage[i]+0.2)) for i, txt in enumerate(labels)]

    ax4.plot(estimated_genome_coverage, shasta_assembled_len, "bo")
    ax4.set(ylabel="Haploid Shasta length", xlabel="estimated read genome coverage",
            title="Estimated Read Genome Cov by Haploid Shasta length")
    [ax4.annotate(txt, (estimated_genome_coverage[i], shasta_assembled_len[i]), xytext=(estimated_genome_coverage[i]+0.2,
                                                shasta_assembled_len[i]+0.2)) for i, txt in enumerate(labels)]
    fig.tight_layout()
    fig.savefig(output_dir + "/shastaLogPlots.png", dpi=150)

def write_out_info(output_dir):
    """ Write all info into log file """
    print("writing data to:", output_dir + "/shasta_stats.csv")
    # if the file is empty write header
    # if not os.path.exists(output_dir + "/shasta_stats.csv"):
    log_out = open(output_dir + "/shasta_stats.csv", 'w')
    fields = ['sample', "read_n50","est_read_cov", "shasta_asm_n50", "shasta_asm_len"]
    log_out.write((',').join(fields) + '\n')
        # log_out.close()

    # open file
    # log_out = open(output_dir + "/shasta_stats.csv", 'a')
    # format line and write
    for i, sam in enumerate(samples):
        outlist = [str(sam), str(read_n50s[i]), str(round(estimated_genome_coverage[i],2)), str(shasta_n50s[i]), str(shasta_assembled_len[i])]
        outstring = (',').join(outlist) + '\n'
        log_out.write(outstring)
    log_out.close()

def main(shastaLog="", shastaLogsFile="", output_dir="output"):

    # if there is a list of files,
    #  loop through that file and run the get_stasta_asm_info function of each file
    # shasta.log.tar.gz

    if len(shastaLogsFile)>0:
        print('extracting files listed in:', shastaLogsFile)
        with open(shastaLogsFile, 'r') as logFiles:
            for file in logFiles.readlines():
                print('reading', file.strip())
                if file.strip()[-6:]=="tar.gz":
                    tfile = tarfile.open(file.strip(), 'r')
                    tfile.extract("shasta.log", "./tfile")
                    get_stasta_asm_info("./tfile/shasta.log")
                    tfile.close()
                else:
                    print('not a tar.gz file')
                    get_stasta_asm_info(file.strip())
                if len(read_n50s) != len(shasta_n50s):
                    print('Problem: lists DIFF LEN', len(read_n50s), read_n50s[-5:], len(shasta_n50s), shasta_n50s[-5:])
                    break
    # for single file input just run get_shasta_asm_info
    if len(shastaLog)>0:
        tfile = tarfile.open(shastaLog, 'r')
        tfile.extract("shasta.log", "./tfile")
        get_stasta_asm_info("./tfile/shasta.log")
        tfile.close()


    if os.path.exists(output_dir):
        print("WARNING: output directory already exists, previous plot will be overwritten.")
    else:
        os.makedirs(output_dir)

    # plot and write data to file
    plotting(output_dir)
    write_out_info(output_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        required=False,
        type=str,
        default="",
        help="Input shasta.log file"
    )

    parser.add_argument(
        "-f",
        required=False,
        type=str,
        default="",
        help="Input file of shasta.log file paths, one on each line"
    )

    parser.add_argument(
        "-o",
        required=False,
        type=str,
        default="output",
        help="Output directory"
    )

    args = parser.parse_args()

    main(shastaLog=args.i, shastaLogsFile=args.f, output_dir=args.o)