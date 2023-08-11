process multiqc {
  container 'quay.io/biocontainers/multiqc:1.15--pyhdfd78af_0'
  memory '32GB'
  cpus 4
  time '1 h'
  publishDir params.outdir, mode: 'copy'

  input:
    file qc_raw
    file trim
    file qc_trim
    file alignment
    file filtered_alignment
    file coverage_alignment
    
  output:
    path "multiqc"

  script:
    """
    mkdir -p multiqc
    multiqc -o multiqc .
    """ 
}
