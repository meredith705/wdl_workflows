version 1.0

workflow runGFAsePhase {
	
	input {
		File assemblyDetailed
		File patKmerFile
		File matKmerFile
		Int kmsize 
		String dockerImage = "meredith705/gfase:latest" 
	}

    # phase gfa
    call gfasePhase {
        input:
            assemblyDetailedGfa = assemblyDetailed,
            patKmerFa=patKmerFile,
            matKmerFa=matKmerFile,
            kSize = kmsize,
            dockerImage=dockerImage
    }

	output {
		File outputMatAssembly = gfasePhase.outMatAssembly
		File outputPatAssembly = gfasePhase.outPatAssembly
	}

}

task gfasePhase {
    input {
		File assemblyDetailedGfa
		File patKmerFa 
		File matKmerFa
		Int kSize 
		# runtime configurations
		Int memSizeGB = 128
		Int threadCount = 1
		Int diskSizeGB = 64
		String dockerImage = "meredith705/gfase:latest"
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
        -k ~{kSize} \
        -i ~{assemblyDetailedGfa} \
        -p ~{patKmerFa} \
        -m ~{matKmerFa} 

    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    output {
        File outMatAssembly = "maternal.fasta"
        File outPatAssembly = "paternal.fasta"
        File outUnphasedAssembly = "unphased.fasta"
        File outPhaseChains = "phase_chains.csv"
    }
}
