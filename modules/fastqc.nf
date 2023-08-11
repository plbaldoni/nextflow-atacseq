process fastqc {
  container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
  memory '32GB'
  cpus 4
  time '1 h'

  input:
    tuple val(sample_id), path(reads)

  output:
    file "*_fastqc.{zip,html}"

  script:
    """
    fastqc --threads $task.cpus --format fastq ${reads[0]} ${reads[1]}
    """
}
