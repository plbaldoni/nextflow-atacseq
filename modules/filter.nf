process filter {
  // ENCODE's post-alignment filtering (read/mate unmapped, not primary, 
  // failing platform, only properly paired) followed by (1) removal of orphan reads 
  // and pairs mapping to different chromosomes, and (2) mark and removal of duplicates
  
  module 'samtools/1.21'
  module 'picard-tools/3.3.0'
  memory '32GB'
  cpus 4
  time '4 h'
  publishDir params.outdir, mode: 'copy', pattern: 'alignment/*.{bam,bam.bai}'
  tag "$sample_id"

  input:
    tuple val(sample_id), path(outbam), path(outstat)

  output:
    tuple val(sample_id), path(nodupbam), path(nodupbambai), path(nodupbam_stat)

  script:
    tmpbam = "${sample_id}.tmp.bam"
    tmpfixbam = "${sample_id}.tmp.fix.bam"
    filtbam = "${sample_id}.filt.bam"
    filttmpbam = "${sample_id}.filt.tmp.bam"
    filtbam_qc="${sample_id}.filt.dup"
    filtbam_stat="${sample_id}.filt.stat"
    
    def opt = !params.singleEnd ? /-F 1804 -f 2/ : /-F 1804/
    
    if( params.trim ) {
      nodupbam ="alignment/${sample_id}.trimmed.filt.nodup.bam"
      nodupbambai = "alignment/${sample_id}.trimmed.filt.nodup.bam.bai"
      nodupbam_stat = "alignment/${sample_id}.trimmed.filt.nodup.stat"
    }
    else {
      nodupbam = "alignment/${sample_id}.filt.nodup.bam"
      nodupbambai = "alignment/${sample_id}.filt.nodup.bam.bai"
      nodupbam_stat = "alignment/${sample_id}.filt.nodup.stat"
    }
    
    """
    mkdir alignment
    
    # Remove read/mate unmapped, not primary alignment, reads failing platform, and low MAPQ reads
    samtools view ${opt} -q ${params.mapq} -u ${outbam} | samtools sort -n /dev/stdin -o ${tmpbam}
    
    # Remove orphan reads and reads maping to different chromosomes
    samtools fixmate -r ${tmpbam} ${tmpfixbam}
    samtools view ${opt} -u ${tmpfixbam} | samtools sort /dev/stdin | samtools addreplacerg -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o ${filtbam} /dev/stdin
    rm -rf ${tmpbam} ${tmpfixbam}
    
    # Marking duplicates and index
    MarkDuplicates \
    --INPUT ${filtbam} \
    --OUTPUT ${filttmpbam} \
    --METRICS_FILE ${filtbam_qc} \
    --VALIDATION_STRINGENCY LENIENT \
    --ASSUME_SORTED true \
    --REMOVE_DUPLICATES false
    mv ${filttmpbam} ${filtbam}
    samtools index  ${filtbam}
    samtools stat  ${filtbam} > ${filtbam_stat}
    
    # Remove duplicates and index
    samtools view ${opt} -b ${filtbam} > ${nodupbam}
    samtools index ${nodupbam}
    samtools stat ${nodupbam} > ${nodupbam_stat}
    """
}
