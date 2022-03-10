hg02.ilmn250.k31.pe.yak version 1.0

workflow runYakAssemblyStats {
	
	input {
		File assemblyFastaPat
		File assemblyFastaMat
		File yakCountPat
		File yakCountMat
		File yakCountSon
		Int shardLinesPerFile = 256000000
		Int fileExtractionDiskSizeGB = 256
		String dockerImage = "juklucas/hpp_yak:latest"

	}

    # get stats
    call yakAssemblyStats {
        input:
            assemblyFastaPat=assemblyFastaPat,
            assemblyFastaMat=assemblyFastaMat,
            patYak=yakCountPat,
            matYak=yakCountMat,
            sampleYak=yakCountSon,
            dockerImage=dockerImage
    }

	output {
		File outputTarball = yakAssemblyStats.outputTarball
		File outputSummary = yakAssemblyStats.outputSummary
	}

}

task yakAssemblyStats {
    input {
        File assemblyFastaPat
        File assemblyFastaMat
        File patYak    # changed to the input yak dbs
        File matYak    # changed also
        File sampleYak
        String genomeSize = "3.2g"
        String minSequenceLength = "100k"
        # runtime configurations
        Int memSizeGB = 128
        Int threadCount = 32
        Int diskSizeGB = 256
        String dockerImage = "juklucas/hpp_yak:latest"
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

        # name
        PREFIX=$(basename ~{assemblyFastaPat} | sed 's/.gz$//' | sed 's/.fa\(sta\)*$//' | sed 's/[._][pm]at\(ernal\)*//')

        # Computing error rates
        yak trioeval -t ~{threadCount} ~{patYak} ~{matYak} ~{assemblyFastaPat} > $PREFIX.pat.yak.switch-error.txt
        yak trioeval -t ~{threadCount} ~{patYak} ~{matYak} ~{assemblyFastaMat} > $PREFIX.mat.yak.switch-error.txt

        # QV
        yak qv -t ~{threadCount} -p -K ~{genomeSize} -l ~{minSequenceLength} ~{sampleYak} ~{assemblyFastaPat} > $PREFIX.pat.yak.qv.txt
        yak qv -t ~{threadCount} -p -K ~{genomeSize} -l ~{minSequenceLength} ~{sampleYak} ~{assemblyFastaMat} > $PREFIX.mat.yak.qv.txt


        # condense
        SUMMARY=$PREFIX.summary.txt
        echo "# mat qv" >>$SUMMARY
        tail -n4 $PREFIX.mat.yak.qv.txt >>$SUMMARY
        echo "# pat qv" >>$SUMMARY
        tail -n4 $PREFIX.pat.yak.qv.txt >>$SUMMARY
        echo "# mat switch" >>$SUMMARY
        tail -n3 $PREFIX.mat.yak.switch-error.txt >>$SUMMARY
        echo "# pat switch" >>$SUMMARY
        tail -n3 $PREFIX.pat.yak.switch-error.txt >>$SUMMARY

        # tar
        tar czvf $PREFIX.yak-qc.tar.gz *txt
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: 1
    }

    output {
        File outputTarball = glob("*.yak-qc.tar.gz")[0]
        File outputSummary = glob("*.summary.txt")[0]
    }
}
