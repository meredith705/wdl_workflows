version 1.0

workflow runLiger2Liger {
	
	input {
		File fastqFile
        Array[File] fastqFiles
		File referenceFile
		String arguments = ""
        String prefix = ""
		Int threads
		String dockerImage = "meredith705/liger2liger:latest"
	}

    call evalChimeras {
		input:
            ref     = referenceFile,
            fastqs  = fastqFile,
            threads = threads,
            args    = arguments,
            prefix  = prefix,
            dockerImage=dockerImage
    }

	output {
		File outTarball = evalChimeras.outTarball
	}

}

task evalChimeras {
    input {
		Array[File] fastqs
		File ref
		Int threads
        String args
        String prefix
		# runtime configurations
		Int memSizeGB = 128
		Int threadCount = threads
		Int diskSizeGB = 64
		String dockerImage
    }
    command <<<
        set -euxo pipefail

        # https://github.com/rlorigro/Liger2LiGer

        # run liger2liger
        evaluate_chimeras.py \
        --ref ~{ref} \
        --fastq ~{sep(",", fastqs)} ~{args} \
        --n_threads ~{threads}

        # tarball
        tar czvf ~{prefix}.liger2liger.tar.gz *.csv *.txt *paf


    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    output {
        File outTarball = prefix + .liger2liger.tar.gz
    }
}
