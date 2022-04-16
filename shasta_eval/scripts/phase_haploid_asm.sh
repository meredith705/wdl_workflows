#!/bin/bash

# Melissa Meredith 04/15/2022
# usage:
# phase_haploid_asm.sh <haploid1.fasta> <haploid2.fasta> <ref_mat.fasta.gz> <ref_pat.fasta.gz>

# Script to bin haploid assemblies into mat, pat, and unphased based off 
#  minimap2 alignments to the reference

# set -e so any non-zero codes will terminate the script
set -e

# store inputs
HAP1=$1
HAP2=$2
# gziped
REF_MAT=$3
REF_PAT=$4

RAW_NUM_HAP1=$(grep -c ">" $HAP1)
RAW_NUM_HAP2=$(grep -c ">" $HAP2)

RAW_NUM_MAT=$(zgrep -c ">" $REF_MAT)
RAW_NUM_PAT=$(zgrep -c ">" $REF_PAT)
TOTAL_NUM_REF_MATPAT=$(($RAW_NUM_MAT + $RAW_NUM_PAT))

# 1 rename ref contigs with mat pat - then combine into one fa
zcat $REF_PAT | sed "s/>/>pat_/g" - > ref.mat.pat.fa
zcat $REF_MAT | sed "s/>/>mat_/g" - >> ref.mat.pat.fa

# 2 rename asm contigs with hap1, hap2 
sed "s/>/>hap1_/g" $HAP1 > hap1.labeled.fa
sed "s/>/>hap2_/g" $HAP2 > hap2.labeled.fa

NUM_REF_CONTIGS=$(grep -c ">" ref.mat.pat.fa)
echo ref $NUM_REF_CONTIGS

NUM_HAP1_LABELED_CONTIGS=$(grep -c ">" hap1.labeled.fa)
NUM_HAP2_LABELED_CONTIGS=$(grep -c ">" hap2.labeled.fa)
echo hap1 $NUM_HAP1_LABELED_CONTIGS
echo hap2 $NUM_HAP2_LABELED_CONTIGS

# confirm that the number of input contigs matches the labeled contigs
# maybe change these to output an error message
[ $RAW_NUM_HAP1 -eq $NUM_HAP1_LABELED_CONTIGS ]
[ $RAW_NUM_HAP2 -eq $NUM_HAP2_LABELED_CONTIGS ]
[ $TOTAL_NUM_REF_MATPAT -eq $NUM_REF_CONTIGS ]
echo raw inputs and labeled intermediates match

# 3 align labled contigs to labled ref
minimap2 \
-a \
-x asm5 \
-k 17 \
-t 64 \
-K 10g \
-I 8g \
--eqx \
--secondary=no \
ref.mat.pat.fa \
hap1.labeled.fa hap2.labeled.fa \
 | samtools view -bh - | samtools sort - -o ref_aligned.hap1.hap2.sorted.bam

# 4 isolate the mat and pat aligned contigs
samtools view ref_aligned.hap1.hap2.sorted.bam | awk '$3~"^mat"' - | awk '$2 == 0 || $2 == 16 {print $1}' - | sort -k1,1 > mat_aligned_contigs.txt
samtools view ref_aligned.hap1.hap2.sorted.bam | awk '$3~"^pat"' - | awk '$2 == 0 || $2 == 16 {print $1}' - | sort -k1,1 > pat_aligned_contigs.txt

# store number of contigs per hap
NUM_MAT_ALIGNED_CONTIGS=$(wc -l mat_aligned_contigs.txt)
NUM_PAT_ALIGNED_CONTIGS=$(wc -l pat_aligned_contigs.txt)

# check for any mat/pat overlap, terminate if any overlaps
NUM_MATPAT_OVERLAP=$( grep -w -f mat_aligned_contigs.txt pat_aligned_contigs.txt | wc -l - )
echo check for overlapping mat/pat contigs
EXP_NUM_OVERLAPS=0
[ $NUM_MATPAT_OVERLAP -eq $EXP_NUM_OVERLAPS ]
echo no overlapping contigs

#  5 identify unphased contigs
echo 
echo identifying unphased contigs
# combine mat pat aligned assembly contigs
cat mat_aligned_contigs.txt > phased.input.contigs.txt
cat pat_aligned_contigs.txt >> phased.input.contigs.txt
# combine all assembly contigs
awk '$0~"^>" {print substr($1,2,length($1))}' hap1.labeled.fa > hap1.2.contigs.txt
awk '$0~"^>" {print substr($1,2,length($1))}' hap2.labeled.fa >> hap1.2.contigs.txt
# isolate contigs not in the phased buckets
grep -v -w -f phased.input.contigs.txt hap1.2.contigs.txt > unphased.aligned_contigs.txt

# 6 extract sequences from hap1, hap2 fastas to make mat, pat fastas
echo
echo extract sequences to make mat pat fastas
# pat
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' pat_aligned_contigs.txt hap1.labeled.fa > pat.haplotype.fa
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' pat_aligned_contigs.txt hap2.labeled.fa >> pat.haplotype.fa
# mat
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' mat_aligned_contigs.txt hap1.labeled.fa > mat.haplotype.fa
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' mat_aligned_contigs.txt hap2.labeled.fa >> mat.haplotype.fa
# unphased
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' unphased.aligned_contigs.txt hap1.labeled.fa > unphased.haplotype.fa
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' unphased.aligned_contigs.txt hap2.labeled.fa >> unphased.haplotype.fa
echo fastas made

echo
echo convert sam to paf for indel analysis
# convert sam to paf for indel analysis
samtools view -h ref_aligned.hap1.hap2.sorted.bam | paftools.js sam2paf - > ref_aligned.hap1.hap2.sorted.paf
