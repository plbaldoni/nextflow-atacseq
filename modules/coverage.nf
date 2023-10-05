process coverage {
  container 'quay.io/biocontainers/deeptools:3.5.1--pyhdfd78af_1'
  memory '64GB'
  cpus 8
  time '1 h'
  publishDir params.outdir, mode: 'copy'
  tag "$sample_id"
  
  input:
    tuple val(sample_id), path(nodupbam), path(nodupbambai), path(nodupbam_stat)

  output:
    tuple val(sample_id), path(outbw)

  script:
  outbw = "coverage/${nodupbam.baseName}.bw"
    """
    mkdir coverage
    bamCoverage --verbose --ignoreForNormalization chrX chrY chrM --normalizeUsing RPGC --binSize 10 \
    --bam ${nodupbam} \
    --outFileName $outbw \
    --outFileFormat bigwig \
    --numberOfProcessors ${task.cpus} \
    --effectiveGenomeSize ${params.gsize}
    """ 
}
