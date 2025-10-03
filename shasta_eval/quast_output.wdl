version 1.0

workflow runQuast {
	input {
        File assemblyFasta
        File? assemblyFasta2
        File? referenceFasta
        File? featuresFile
        String extraArguments="--large --no-html --no-icarus --no-plots"
        Int memSizeGB = 64
        Int threadCount = 16
        Int diskSizeGB = 64
        String dockerImage = "meredith705/quast-5.3.0:latest"
    }

    call quast {
        input:
            assemblyFasta = assemblyFasta,
            assemblyFasta2 = assemblyFasta2,
            referenceFasta = referenceFasta,
            extraArguments = extraArguments,
            featuresFile = featuresFile,
            memSizeGB = memSizeGB,
            threadCount = threadCount,
            diskSizeGB = diskSizeGB,
            dockerImage = dockerImage
    }

    output {
        File outTarball = quast.outputTarball
        File outSummary = quast.outputSummary

    }
}

task quast {
    input {
        File assemblyFasta
        File? assemblyFasta2
        File? referenceFasta
        File? featuresFile
        String extraArguments="--large --no-html --no-icarus --no-plots"
        Int memSizeGB = 64
        Int threadCount = 16
        Int diskSizeGB = 64
        String dockerImage = "tpesout/hpp_quast:latest"
    }

	command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        # initilization
        ASM_FILENAME=$(basename -- "~{assemblyFasta}")
        if [[ $ASM_FILENAME =~ \.gz$ ]]; then
            cp ~{assemblyFasta} .
            gunzip $ASM_FILENAME
            ASM_FILENAME="${ASM_FILENAME%.gz}"
        else
            ln -s ~{assemblyFasta}
        fi
        PREFIX="${ASM_FILENAME%.*}"

        # init quast command
        cmd=(python /opt/quast/quast-5.3.0/quast-lg.py )
        cmd+=( -t ~{threadCount} )
        cmd+=( -o $PREFIX.quast )

        
        # include reference fasta if supplied
        if [[ -f "~{referenceFasta}" ]]; then
            REF_FILENAME=$(basename -- "~{referenceFasta}")
            if [[ $REF_FILENAME =~ \.gz$ ]]; then
                cp ~{referenceFasta} .
                gunzip $REF_FILENAME
                REF_FILENAME="${REF_FILENAME%.gz}"
            else
                ln -s ~{referenceFasta}
            fi
            cmd+=( -r $REF_FILENAME )
        fi

        # include features ( like genes ) if provided
        if [[ -f "~{featuresFile}" ]]; then
            FEATURE_FILENAME=$(basename -- "~{featuresFile}")
            ln -s ~{featuresFile}
            cmd+=( --features gene:$FEATURE_FILENAME )  ### do we want just genes? 
        fi


        # include extra arguments if supplied
        if [[ ! -z "~{extraArguments}" ]]; then
            cmd+=( ~{extraArguments} )
        fi

        # finalize command with first assembly haplotype
        cmd+=( $ASM_FILENAME )

        # include second haplotype assembly fasta if supplied
        if [[ -f "~{assemblyFasta2}" ]]; then
            ASM2_FILENAME=$(basename -- "~{assemblyFasta2}")
            if [[ $ASM2_FILENAME =~ \.gz$ ]]; then
                cp ~{assemblyFasta2} .
                gunzip $ASM2_FILENAME
                ASM2_FILENAME="${ASM2_FILENAME%.gz}"
            else
                ln -s ~{assemblyFasta2}
            fi
            cmd+=( $ASM2_FILENAME )
        fi

        # run command
        "${cmd[@]}"

        # save output
        tar czvf $PREFIX.quast.tar.gz $PREFIX.quast

	>>>
	output {
		File outputTarball = glob("*.quast.tar.gz")[0]
		File outputSummary = glob("*.quast/report.txt")[0]
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
    }
}
