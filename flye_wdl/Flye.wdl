version 1.0


workflow Flye {

    meta {
        description: "Assemble a genome using Flye"
    }
    parameter_meta {
        genome_size: "Estimated genome size in base pairs"
        reads: "Input reads (in fasta or fastq format, compressed or uncompressed)"
        prefix: "Prefix to apply to assembly output filenames"
    }

    input {
        File reads
        String prefix
    }

    call Assemble {
        input:
            reads  = reads,
            prefix = prefix,
    }

    output {
        File gfa = Assemble.gfa
        File fa = Assemble.fa
    }
}

task Assemble {
    input {
        File reads
        String prefix = "out"
        String flye_args = "--nano-hq"
        String dockerImage = "quay.io/biocontainers/flye:2.9.1--py310h590eda1_0"
        Int threads = 24
        Int memSizeGB = 768
        Int disk_size = 1125
        Int preemptible_tries = 0
        Int max_retries = 0
    }

    parameter_meta {
        reads:    "reads (in fasta or fastq format, compressed or uncompressed)"
        prefix:   "prefix to apply to assembly output filenames"
    }


    command <<<
        set -euxo pipefail


        flye ~{flye_args} ~{reads} --threads ~{threads} --out-dir asm

        mv asm/assembly.fasta ~{prefix}.flye.fa
        mv asm/assembly_graph.gfa ~{prefix}.flye.gfa
    >>>

    output {
        File fa = "~{prefix}.flye.fa"
        File gfa = "~{prefix}.flye.gfa"
    }


    runtime {
        cpu:                    threads
        memory:                 memSizeGB + " GB"
        disks: "local-disk " +  disk_size + " LOCAL"
        preemptible:            preemptible_tries
        maxRetries:             max_retries
        docker:                 dockerImage
    }
}
