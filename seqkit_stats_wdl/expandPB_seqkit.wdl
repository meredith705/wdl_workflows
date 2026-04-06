version 1.0

workflow expandTar {

  meta {
    description: "expand PacBio Tars."
  }

  input {
    File   tar_file
    String docker = "meredith705/seqkitstats:latest"
  }

  call expand {
    input:
      tar_file       = tar_file,
      docker         = docker
  }

  call seqkitStats{
    input:
      fastq          = expand.hifi_fastq
  }

  output {
    File  hifi_read_bam      = expand.hifi_read_bam
    File  hifi_read_bam_pbi  = expand.hifi_read_bam_pbi
    File  seqkitStatsReport  = seqkitStats.outStats
  }
}


task expand {
  input {
    File tar_file
    String docker
    Int cpu = 2
    Int memory_gb = 96
    Int disk_gb = 150
    Int preemptible = 1
  }

  String tar_basename = sub(basename(tar_file), "\\_cell1.tar", "")

  command <<<
    set -o pipefail
    set -e
    set -u
    set -o xtrace

    # expand the tar
    tar -xvf ~{tar_file}

    ls

    echo "~{tar_basename}"

  >>>

  output {
    File  hifi_read_bam      = "~{tar_basename}.hifi_reads.bam"
    File  hifi_read_bam_pbi  = "~{tar_basename}.hifi_reads.bam.pbi"
    File  hifi_fastq      = "~{tar_basename}_HiFi.fastq.gz"
  }

  runtime {
    docker:  docker
    cpu:     "~{cpu}"
    memory:  "~{memory_gb} GB"
    disks:   "local-disk ~{disk_gb} HDD"
    preemptible: "~{preemptible}"
  }
}



task seqkitStats {
    input {
      File fastq
      # runtime configurations
      Int memSizeGB = 128
      Int threadCount = 10
      Int diskSizeGB = 96
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

        seqkit stats -j 10 -a ~{fastq} -T > ~{fastq}.stats.tsv

    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    output {
        File outStats = "~{fastq}.stats.tsv"

    }
}