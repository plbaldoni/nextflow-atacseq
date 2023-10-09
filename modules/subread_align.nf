process subread_align {
  module 'subread/2.0.6'
  module 'samtools/1.18'
  memory '64GB'
  cpus 8
  time '8 h'
  tag "$sample_id"

  input:
    tuple val(sample_id), path(reads)
    
  output:
    tuple val(sample_id), path(outbam), path(outstat)

  script:
    outbam = "${sample_id}.bam"
    outstat = "${sample_id}.stat"
    def single = reads instanceof Path
    def read1 = !single ? /-r "${reads[0]}"/ : /-r "${reads}"/
    def read2 = !single ? /-R "${reads[1]}"/ : ''
    """
    subread-align --sortReadsByCoordinates --multiMapping \
    -T $task.cpus \
    -i $params.subreadIndex \
    -t 1 \
    -D $params.maxfraglen \
    -B $params.maxmultimap \
    ${read1} ${read2} \
    -o ${outbam}
    
    samtools index ${outbam}
    samtools stat ${outbam} > ${outstat}
    """
}
