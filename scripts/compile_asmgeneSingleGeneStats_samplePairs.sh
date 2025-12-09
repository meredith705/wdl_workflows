#!/usr/bin/env bash

# directories containing shasta and hifiasm Asmgene files
dir_s="/Users/mmmeredi/Desktop/paten/ROSMAP_brain_samples/asm_quast/asmgene/AsmgeneStats/perGeneStats/AsmgenePerGeneStats"
dir_h="/Users/mmmeredi/Desktop/paten/ROSMAP_brain_samples/asm_quast/asmgene/hifiasmAsmgeneStats/AsmgenePerGeneStats/hifiasmAsmgenePerGeneStats"

# loop through all Shasta files
for sfile in "$dir_s"/*per_gene_stats.txt; do
    # extract the sample base name (e.g., sampleB from sampleB_s.txt)
    sample=$(basename "$sfile")
    sample=${sample%%_hapdup_dual_*}

    # find corresponding _hifiasm. file
    hfile="$dir_h/${sample}_hifiasm.ont.bp.hap1.p_ctg.0.97.per_gene_stats.txt"

    # check if both files exist
    if [[ -f "$sfile" && -f $hfile ]]; then
        echo "python script.py -a \"$sfile\" -b \"$hfile\" -o \"asmgenePerStats_${sample}.tsv\""
        # uncomment to actually run:
        python3 ~/Desktop/paten/wdl_git_workflows/wdl_workflows/scripts/compareAsmgenePerGeneStats.py -a "$sfile" -b "$hfile" -o "asmgenePerStats_${sample}.tsv"
    else
        echo "Warning: missing file for sample $sample" >&2
        echo "shasta file: $sfile"
        echo "hifiasm file: $hfile"
    fi
done
