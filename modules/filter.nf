process filter {
  // ENCODE's post-alignment filtering (read/mate unmapped, not primary, 
  // failing platform, only properly paired) followed by (1) removal of orphan reads 
  // and pairs mapping to different chromosomes, and (2) mark and removal of duplicates
  
  module 'samtools/1.17'
  module 'picard-tools/2.26.11'
  memory '32GB'
  cpus 4
  time '4 h'
  publishDir params.outdir, mode: 'copy', pattern: '*.{filt.nodup.bam,filt.nodup.bam.bai}'


  input:
    file alignment

  output:
    file "*.{filt.nodup.bam,filt.nodup.bam.bai,filt.nodup.stat,.filt.dup,.filt.stat}"

  script:
    bam = "${alignment[0]}"
    tmpbam = "${alignment[0].baseName}.tmp.bam"
    tmpfixbam = "${alignment[0].baseName}.tmp.fix.bam"
    filtbam = "${alignment[0].baseName}.filt.bam"
    filttmpbam = "${alignment[0].baseName}.filt.tmp.bam"
    filtbam_qc="${alignment[0].baseName}.filt.dup"
    filtbam_stat="${alignment[0].baseName}.filt.stat"
    nodupbam="${alignment[0].baseName}.filt.nodup.bam"
    nodupbam_stat="${alignment[0].baseName}.filt.nodup.stat"
    
    """
    samtools view -F 1804 -f 2 -q ${params.mapq} -u ${bam} | samtools sort -n /dev/stdin -o ${tmpbam}
    samtools fixmate -r ${tmpbam} ${tmpfixbam}
    samtools view -F 1804 -f 2 -u ${tmpfixbam} | samtools sort /dev/stdin -o ${filtbam}
    rm -rf ${tmpbam} ${tmpfixbam}
    
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
    
    samtools view -F 1804 -f 2 -b ${filtbam} > ${nodupbam}
    samtools index ${nodupbam}
    samtools stat ${nodupbam} > ${nodupbam_stat}
    """
}
