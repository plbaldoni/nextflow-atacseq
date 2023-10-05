process trimgalore {
  container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'
  memory '32GB'
  cpus 4
  time '1 h'
  tag "$sample_id"

  input:
    tuple val(sample_id), path(reads)
    
  output:
    tuple val(sample_id), path("*.gz")
    file "*.txt"
    
  script:
    def single = reads instanceof Path
    def read1 = !single ? /"--paired ${reads[0]}"/ : /"${reads}"/
    def read2 = !single ? /"${reads[1]}"/ : ''
    
    """
    trim_galore -e $params.erate --cores $task.cpus ${read1} ${read2}
    """
}
