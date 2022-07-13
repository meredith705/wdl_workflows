# seqkit stats 

version 1.0

workflow runSeqkitStats {
	
	input {
		File inputPassFastq
		File inputFailFastq
		String dockerImage = "meredith705/seqkitstats:latest" 
	}

    call seqkitStats {
        input:
            passFastq = inputPassFastq,
            failFastq = inputFailFastq,
            dockerImage=dockerImage
    }

	output {
		File outputPassStats = seqkitStats.outPassStats
		File outputFailStats = seqkitStats.outFailStats
	}

}

task seqkitStats {
    input {
		File passFastq
		File failFastq 
		# runtime configurations
		Int memSizeGB = 128
		Int threadCount = 10
		Int diskSizeGB = 64
		String dockerImage = "meredith705/seqkitstats:latest"
    }
    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        seqkit stats -j 10 -a ~{passFastq} -T > ~{passFastq}.stats.pass.tsv
        seqkit stats -j 10 -a ~{failFastq} -T > ~{failFastq}.stats.fail.tsv 

    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    output {
        File outPassStats = "~{passFastq}.stats.pass.tsv"
        File outFailStats = "~{failFastq}.stats.fail.tsv"

    }
}