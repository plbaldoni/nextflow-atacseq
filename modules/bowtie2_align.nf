process bowtie2_align {
  container 'quay.io/biocontainers/bowtie2:2.5.1--py38he00c5e5_2'
  container 'quay.io/biocontainers/samtools:1.17--hd87286a_1'
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
