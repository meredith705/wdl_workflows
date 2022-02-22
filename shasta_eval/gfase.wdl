version 1.0

workflow runGFAsePhase {
	
	input {
		File AssemblyDetailed
        File patKmerFile
        File matKmerFile
        String dockerImage = "meredith705/gfase:latest" # parsing error here

	}

    # phase gfa
    call gfasePhase {
        input:
            assemblyDetailedGfa = AssemblyDetailed,
            patKmerFa=patKmerFile,
            matKmerFa=matKmerFile,
            dockerImage=dockerImage
    }

	output {
		#File outputTarball = yakAssemblyStats.outputTarball
		#File outputSummary = yakAssemblyStats.outputSummary
        File outputMatAssembly = gfasePhase.matAssembly
        File outputPatAssembly = gfasePhase.patAssembly
	}

}

task gfasePhase {
    input {
        File assemblyDetailedGfa
        File patKmerFa 
        File matKmerFa
        Int kSize = 31
        # runtime configurations
        Int memSizeGB = 128
        Int threadCount = 1
        Int diskSizeGB = 64
        String dockerImage = ""meredith705/gfase:latest""
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

        # Phase the gfa
        phase_haplotype_paths \
        -i {assemblyDetailedGfa} \
        -p {patKmerFa} \
        -m {matKmerFa} \
        -k {kSize}

        # name
        #PREFIX=$(basename ~{assemblyFastaPat} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')

    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: 1
    }

    output {
        File outputMatAssembly = maternal.fasta
        File outputPatAssembly = paternal.fasta
    }
}
