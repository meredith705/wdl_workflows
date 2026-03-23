version 1.0

# counts the total number of reads and the number of unmapped reads
# in a BAM file using samtools flagstat and samtools view.


workflow CountReads {

  meta {
    description: "Count total and unmapped reads in a BAM file using samtools."
    author: "Melissa Meredith, Datatecnica LLC"
  }

  input {
    File   bam_file
    File   bam_index          
    String samtools_docker = "biocontainers/samtools:v1.9-4-deb_cv1"
  }

  call FlagstatAndCount {
    input:
      bam_file       = bam_file,
      bam_index      = bam_index,
      docker         = samtools_docker,

  }

  output {
    Int  total_reads     = FlagstatAndCount.total_reads
    Int  primary_reads   = FlagstatAndCount.primary_reads
    Int  unmapped_reads  = FlagstatAndCount.unmapped_reads
    File flagstat_report = FlagstatAndCount.flagstat_report
  }
}


task FlagstatAndCount {

  meta {
    description: "Run samtools flagstat and parse total / unmapped read counts."
  }

  input {
    File   bam_file
    File   bam_index
    String docker
    Int    cpu       = 16
    Int    memory_gb = 96
    Int    disk_gb   = round(size(bam_file, 'G')) + 30
    Int preemptible  = 0
  }

  # Localise the index next to the BAM when supplied
  String bam_basename = basename(bam_file)

  command <<<
    set -euo pipefail

    # samtools flagstat 
    samtools flagstat -@~{cpu} ~{bam_file} > flagstat.txt
    echo "=== flagstat output ===" >&2
    cat flagstat.txt >&2

    # Parse total reads from flagstat output
    # "<N> + <N> in total (QC-passed reads + QC-failed reads)"
    TOTAL=$(grep -m1 "in total" flagstat.txt | awk '{print $1 + $3}')
    PRIMARY=$(grep -m1 "primary" flagstat.txt | awk '{print $1 + $3}')

    # Count unmapped reads via samtools view -c -f 4 
    UNMAPPED=$(samtools view -@~{cpu} -c -f 4 ~{bam_file})

    echo "$TOTAL"    > total_reads.txt
    echo "$PRIMARY"  > primary_reads.txt
    echo "$UNMAPPED" > unmapped_reads.txt

    echo "Total reads   : $TOTAL"    >&2
    echo "Unmapped reads: $UNMAPPED" >&2
  >>>

  output {
    Int  total_reads     = read_int("total_reads.txt")
    Int  primary_reads   = read_int("primary_reads.txt")
    Int  unmapped_reads  = read_int("unmapped_reads.txt")
    File flagstat_report = "flagstat.txt"
  }

  runtime {
    docker:  docker
    cpu:     "~{cpu}"
    memory:  "~{memory_gb} GB"
    disks:   "local-disk ~{disk_gb} HDD"
    preemptible: "~{preemptible}"
  }
}