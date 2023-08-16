process bowtie2_align {
  module 'bowtie2/2.4.4'
  module 'samtools/1.18'
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
    bowtie2 --mm \
    -k $params.maxmultimap \
    -X $params.maxfraglen \
    --threads $task.cpus \
    -x $params.bowtie2Index \
   	-1 ${reads[0]} -2 ${reads[1]} |
   	samtools view -Su /dev/stdin | samtools sort - -o ${outbam}
    
    samtools index ${outbam}
    samtools stat ${outbam} > ${outstat}
    """
}
