process bowtie2_align {
  module 'bowtie2/2.5.3'
  module 'samtools/1.21'
  memory '72GB'
  cpus 12
  time params.bowtie2Time
  tag "$sample_id"

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path(outbam), path(outstat)

  script:
    outbam = "${sample_id}.bam"
    outstat = "${sample_id}.stat"
    def single = reads instanceof Path
    def read1 = !single ? /-1 "${reads[0]}"/ : /-U "${reads}"/
    def read2 = !single ? /-2 "${reads[1]}"/ : ''
    """
    bowtie2 --mm ${params.optBowtie2} \
    -x ${params.bowtie2Index} \
    --threads ${task.cpus} \
    ${read1} ${read2} |
    samtools view -Su /dev/stdin | samtools sort - -o ${outbam}

    samtools index ${outbam}
    samtools stat ${outbam} > ${outstat}
    """
}
