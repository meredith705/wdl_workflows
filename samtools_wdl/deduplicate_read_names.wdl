version 1.0

workflow DeduplicateReadNames {
    meta {
        description: "Appends a unique line number to each read name in a BAM file to resolve duplicate read names after merging multiple BAMs."
    }

    input {
        File input_bam
        String sample_name
        String docker_image = "meredith705/card_sniffles:2.7.2"
    }

    call DeduplicateReadNamesTask {
        input:
            input_bam    = input_bam,
            sample_name  = sample_name,
            docker_image = docker_image
    }

    output {
        File output_bam     = DeduplicateReadNamesTask.output_bam
        File output_bam_bai = DeduplicateReadNamesTask.output_bam_bai
    }
}

task DeduplicateReadNamesTask {
    meta {
        description: "Uses samtools + awk to append the SAM record line number to each read name, guaranteeing uniqueness."
    }

    input {
        File   input_bam
        String sample_name
        String docker_image

        # Runtime
        Int    cpu          = 4
        Int    mem_gb       = 16
        Int    disk_gb      = 3 * ceil(size(input_bam, "GB")) + 20
        Int    preemptible  = 1
    }

    command <<<
        set -euo pipefail

        echo "INFO: Input BAM: ~{input_bam}"
        echo "INFO: Sample: ~{sample_name}"

        # Stream BAM -> SAM, rename reads, convert back to sorted BAM
        samtools view -h ~{input_bam} \
            | awk 'BEGIN { OFS="\t" }
                   /^@/ { print; next }
                   { $1 = $1 "_" NR; print }' \
            | samtools sort \
                -@ ~{cpu} \
                -m 2G \
                -o ~{sample_name}.dedup_readnames.bam

        # Index the output
        samtools index ~{sample_name}.dedup_readnames.bam

        # Quick sanity check — verify no duplicate read names remain
        echo "INFO: Checking for remaining duplicates (should be 0)..."
        DUPES=$(samtools view ~{sample_name}.dedup_readnames.bam \
                    | awk '{print $1}' \
                    | sort | uniq -d | wc -l)
        echo "INFO: Duplicate read names remaining: $DUPES"
        if [ "$DUPES" -gt 0 ]; then
            echo "ERROR: $DUPES duplicate read names still present!" >&2
            exit 1
        fi

        echo "INFO: Done."
    >>>

    output {
        File output_bam     = "~{sample_name}.dedup_readnames.bam"
        File output_bam_bai = "~{sample_name}.dedup_readnames.bam.bai"
    }

    runtime {
        docker:      docker_image
        cpu:         cpu
        memory:      mem_gb + " GB"
        disks:       "local-disk " + disk_gb + " SSD"
        preemptible: preemptible
    }
}
