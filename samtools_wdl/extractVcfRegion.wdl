version 1.0

workflow ExtractRegionFromGVCFs {
    input {
        Array[File] gvcf_files
        Array[File] gvcf_indices
        String region                  # ex: "chr15:34426000-34427000"
        String output_prefix           
        File reference_fasta
        File reference_fasta_fai
    }

    # Extract region from each gVCF individually
    scatter (i in range(length(gvcf_files))) {
        call ExtractRegion {
            input:
                gvcf        = gvcf_files[i],
                gvcf_index  = gvcf_indices[i],
                region      = region,
        }
    }

    # Merge all per-sample region VCFs into one
    call MergeAndCount {
        input:
            vcf_files     = ExtractRegion.region_vcf,
            vcf_indices   = ExtractRegion.region_vcf_index,
            output_prefix = output_prefix,
            region        = region,
            reference_fasta     = reference_fasta,
            reference_fasta_fai = reference_fasta_fai,
    }

    output {
        File merged_vcf       = MergeAndCount.merged_vcf
        File merged_vcf_index = MergeAndCount.merged_vcf_index
        File variant_counts   = MergeAndCount.variant_counts
    }
}

task ExtractRegion {
    input {
        File   gvcf
        File   gvcf_index
        String region

        Int    disk_gb   = 50
        Int    memory_gb = 4
        Int    cpu       = 1
        String docker    = "meredith705/truvari"
    }

    String sample_name = basename(gvcf, ".g.vcf.gz")

    command <
        set -euo pipefail

        bcftools view \
            -r ~{region} \
            -O z \
            -o ~{sample_name}.~{region}.vcf.gz \
            ~{gvcf}

        bcftools index -t ~{sample_name}.~{region}.vcf.gz
    >>>

    output {
        File region_vcf       = "~{sample_name}.~{region}.vcf.gz"
        File region_vcf_index = "~{sample_name}.~{region}.vcf.gz.tbi"
    }

    runtime {
        docker:   docker
        cpu:      cpu
        memory:   "~{memory_gb} GB"
        disks:    "local-disk ~{disk_gb} HDD"
    }
}

task MergeAndCount {
    input {
        Array[File] vcf_files
        Array[File] vcf_indices
        String      output_prefix
        String      region
        File        reference_fasta
        File        reference_fasta_fai

        Int    disk_gb   = 100
        Int    memory_gb = 8
        Int    cpu       = 2
        String docker    = "meredith705/truvari"
    }

    command <
        set -euo pipefail

        # Write VCF list to file
        VCF_LIST="vcf_list.txt"
        echo "~{sep='\n' vcf_files}" > $VCF_LIST

        # Merge all per-sample VCFs
        bcftools merge \
            --file-list $VCF_LIST \
            --regions ~{region} \
            --output-type z \
            --output ~{output_prefix}.vcf.gz

        bcftools index -t ~{output_prefix}.vcf.gz

        # Count variants by type
        OUTFILE="~{output_prefix}.variant_counts.txt"

        echo "Region: ~{region}" > $OUTFILE
        echo "Output: ~{output_prefix}.vcf.gz" >> $OUTFILE
        echo "" >> $OUTFILE

        TOTAL=$(bcftools view -H ~{output_prefix}.vcf.gz | wc -l)
        SNPS=$(bcftools view -H -v snps ~{output_prefix}.vcf.gz | wc -l)
        INDELS=$(bcftools view -H -v indels ~{output_prefix}.vcf.gz | wc -l)
        MNPS=$(bcftools view -H -v mnps ~{output_prefix}.vcf.gz | wc -l)
        OTHER=$(bcftools view -H -v other ~{output_prefix}.vcf.gz | wc -l)

        echo "Variant counts" >> $OUTFILE
        echo "--------------" >> $OUTFILE
        printf "%-10s %s\n" "Total:"  "$TOTAL"  >> $OUTFILE
        printf "%-10s %s\n" "SNPs:"   "$SNPS"   >> $OUTFILE
        printf "%-10s %s\n" "Indels:" "$INDELS" >> $OUTFILE
        printf "%-10s %s\n" "MNPs:"   "$MNPS"   >> $OUTFILE
        printf "%-10s %s\n" "Other:"  "$OTHER"  >> $OUTFILE

        echo "" >> $OUTFILE

        # Per-sample counts
        echo "Per-sample non-ref genotype counts" >> $OUTFILE
        echo "------------------------------------" >> $OUTFILE
        bcftools stats -s - ~{output_prefix}.vcf.gz \
            | grep "^PSC" \
            | awk 'BEGIN{printf "%-30s %-10s %-10s\n", "Sample", "Hom_RR", "Het"} \
                   {printf "%-30s %-10s %-10s\n", $3, $5, $6}' \
            >> $OUTFILE
    >>>

    output {
        File merged_vcf       = "~{output_prefix}.vcf.gz"
        File merged_vcf_index = "~{output_prefix}.vcf.gz.tbi"
        File variant_counts   = "~{output_prefix}.variant_counts.txt"
    }

    runtime {
        docker:   docker
        cpu:      cpu
        memory:   "~{memory_gb} GB"
        disks:    "local-disk ~{disk_gb} HDD"
    }
}