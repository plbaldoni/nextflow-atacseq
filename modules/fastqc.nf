process fastqc {
  container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
  memory '32GB'
  cpus 4
  time '1 h'
  tag "$sample_id"

  input:
    tuple val(sample_id), path(reads)

  output:
    file "*_fastqc.{zip,html}"

  script:
    def single = reads instanceof Path
    def read1 = !single ? /"${reads[0]}"/ : /"${reads}"/
    def read2 = !single ? /"${reads[1]}"/ : ''
    """
    fastqc --threads $task.cpus --format fastq ${read1} ${read2}
    """
}
