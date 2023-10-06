process trimgalore {
  container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'
  memory '32GB'
  cpus 4
  time '1 h'
  tag "$sample_id"

  input:
    tuple val(sample_id), path(reads)
    
  output:
    tuple val(sample_id), path("*.trimmed.fastq.gz")
    file "*.txt"
    
  script:
    def single = reads instanceof Path
    def read1 = !single ? /--paired "${reads[0]}"/ : /"${reads}"/
    def read2 = !single ? /"${reads[1]}"/ : ''
    
    if( params.singleEnd ) {
      mvr1 = /mv *val*.fq.gz "${reads.simpleName}".trimmed.fastq.gz/
      mvr2 = ''
    }
    else {
      mvr1 = /mv *_val_1.fq.gz "${reads[0].simpleName}".trimmed.fastq.gz/
      mvr2 = /mv *_val_2.fq.gz "${reads[1].simpleName}".trimmed.fastq.gz/
    }
    """
    trim_galore -e $params.erate --cores $task.cpus ${read1} ${read2}
    ${mvr1}
    ${mvr2}
    """
}
