#!/bin/bash
#SBATCH --job-name=hifiasm
#SBATCH --mail-type=FAIL,END
#SBATCH --partition=largemem
#SBATCH --mail-user=melissa@datatecnica.com
#SBATCH --nodes=1
#SBATCH --mem=1000G
#SBATCH --cpus-per-task=96
#SBATCH --output=%x.%j.log
#SBATCH --time=12:00:00

# Usage: sbatch --mail-type=BEGIN,TIME_LIMIT_90,END script.sh sample_name file1.bam file2.bam file3.bam

set -eox

# go to directory 
cd ...
date
hostname
pwd

module load hifiasm
module load samtools 

sample_name=$1
shift

for READS in "$@"; do
    samtools fastq -@ 96 $READS >> ${sample_name}.merged.fastq
done



# run hifiasm
hifiasm -t $SLURM_CPUS_PER_TASK -o $sample_name.hifiasm.ont ${sample_name}.merged.fastq 2> $sample_name.hifiasm.log


awk '/^S/{print ">"$2;print $3}' $sample_name.hifiasm.ont.bp.hap1.p_ctg.gfa > $sample_name.hifiasm.ont.bp.hap1.p_ctg.fa
awk '/^S/{print ">"$2;print $3}' $sample_name.hifiasm.ont.bp.hap2.p_ctg.gfa > $sample_name.hifiasm.ont.bp.hap2.p_ctg.fa

bgzip -@ 64 $sample_name.hifiasm.ont.bp.hap1.p_ctg.fa
bgzip -@ 64 $sample_name.hifiasm.ont.bp.hap2.p_ctg.fa

bgzip -@ 96 ${sample_name}.merged.fastq

