process coverage {
  container 'quay.io/biocontainers/deeptools:3.5.1--pyhdfd78af_1'
  memory '64GB'
  cpus 8
  time '1 h'
  publishDir params.outdir, mode: 'copy'
  
  input:
    file alignment

  output:
    file "${alignment[0].baseName}.bw"

  script:
    """
    bamCoverage \
    --bam ${alignment[0]} \
    --outFileName ${alignment[0].baseName}.bw \
    --outFileFormat bigwig \
    --numberOfProcessors ${task.cpus} \
    --effectiveGenomeSize ${params.gsize} \
    --normalizeUsing CPM \
    --ignoreForNormalization chrX chrM \
    --verbose
    """ 
}
