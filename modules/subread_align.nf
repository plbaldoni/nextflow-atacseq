process subread_align {
  module 'subread/2.0.6'
  module 'samtools/1.17'
  memory '64GB'
  cpus 8
  time '4 h'

  input:
    tuple val(sample_id), path(reads)
    
  output:
    file "*.{bam,bai,stat}"

  script:
    outbam = "${sample_id}.bam"
    outstat = "${sample_id}.stat"
    """
    subread-align --sortReadsByCoordinates --multiMapping \
    -T $task.cpus \
    -i $params.subreadIindex \
    -t 1 \
    -D $params.maxfraglen \
    -B $params.maxmultimap \
    -r ${reads[0]} \
    -R ${reads[1]} \
    -o ${outbam}
    
    samtools index ${outbam}
    samtools stat ${outbam} > ${outstat}
    """
}
