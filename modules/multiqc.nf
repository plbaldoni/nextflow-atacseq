process multiqc {
  module 'MultiQC/1.10.1'
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
