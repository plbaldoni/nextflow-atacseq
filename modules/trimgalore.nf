process trimgalore {
  module 'trimgalore/0.6.7'
  memory '32GB'
  cpus 4
  time '1 h'

  input:
    tuple val(sample_id), path(reads)
    
  output:
    tuple val(sample_id), path("*.gz")
    file "*.txt"
    
  script:
    """
    trim_galore --paired -e $params.erate --cores $task.cpus ${reads[0]} ${reads[1]}
    """
}
