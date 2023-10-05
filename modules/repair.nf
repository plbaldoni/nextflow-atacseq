process repair {
  container 'quay.io/biocontainers/bbmap:39.01--h92535d8_1'
  memory '128GB'
  cpus 4
  time '1 h'
  tag "$sample_id"

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("*.fastq")

  script:
    """
    zcat ${reads[0]} > in1.fq
    zcat ${reads[1]} > in2.fq
    repair.sh -Xmx100g in=in1.fq in2=in2.fq out=${reads[0].simpleName}.fastq out2=${reads[1].simpleName}.fastq
    """
}
