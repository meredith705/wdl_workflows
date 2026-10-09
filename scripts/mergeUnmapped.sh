#!/bin/bash
#SBATCH --job-name=margin_wg2
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=long
#SBATCH --mail-user=mmmeredi@ucsc.edu
#SBATCH --nodes=1
#SBATCH --mem=500G
#SBATCH --cpus-per-task=64
#SBATCH --output=%x.%j.log
#SBATCH --time=72:00:00

set -eox

# go to directory
cd /private/groups/patenlab/memeredith/CARD_PPMI/margin
date;hostname;pwd

unphasedMappedBAM=$1
phasedBAM=$2
sample=$3
threads=64


outname="${sample}"+".GRCh38.bam"


#1 : get the header from the first unphased bam into a tmp.sam to append all the reads to
samtools view -H "${unphasedMappedBAM}" > tmp.extracted_reads.sam

# 2: get unmapped reads from bam
echo "find unmapped reads"
samtools view -f 4 -@ "${threads}" "${unphasedMappedBAM}" >> tmp.extracted_reads.sam


UNMAPPED=$(samtools view -c -f 4 -@ "${threads}" tmp.extracted_reads.sam )
#"${unphasedMappedBAM}")
echo "Unmapped reads: ${UNMAPPED}"
echo "Unmapped reads: $UNMAPPED" >&2
echo "$UNMAPPED" > readcount.txt

# 3: convert the tmp.sam to a bam
samtools view -b -@ "${threads}" tmp.extracted_reads.sam | samtools sort -@ "${threads}" - > tmp.extracted_reads.bam

# 4: merge unmapped and haplotagged bams
samtools merge -@ "${threads}" -o "${outname}" "${phasedBAM}" tmp.extracted_reads.bam

# 5: index the merged BAM
samtools index -@ "${threads}" "${outname}"
