params.reads = "$launchDir/*{R1,R2}_001.fastq.gz"
params.outdir = "$launchDir"
params.erate = 0.2
params.mapq = 30
params.gsize = 2654605538 // GRCm39 genome size (https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001635.27/)
params.maxmultimap = 4
params.maxfraglen = 2000
params.subreadIndex = "$launchDir/subread_index"
params.bowtie2Index = "$launchDir/bowtie2_index"
params.subread = false
params.bowtie2 = false
params.repair = false
params.trim = false
params.singleEnd = false
params.bowtie2Time = "24h"
params.subreadTime = "24h"

include { fastqc } from './modules/fastqc'
include { fastqc as fastqc_trim} from './modules/fastqc'
include { trimgalore } from './modules/trimgalore'
include { bowtie2_align } from './modules/bowtie2_align'
include { subread_align } from './modules/subread_align'
include { filter } from './modules/filter'
include { coverage } from './modules/coverage'
include { multiqc_trim } from './modules/multiqc_trim'
include { multiqc } from './modules/multiqc'
include { repair } from './modules/repair'

workflow {
  
  if ( params.repair ) {
    if ( params.singleEnd ) {
      error "Repair is only possible for paired-end reads"
    }
    else {
      ch_reads_input = repair(Channel.fromFilePairs(params.reads, checkIfExists: true)) 
    }
  }
  else {
    ch_reads_input = Channel.fromFilePairs(params.reads, size: params.singleEnd ? 1 : 2, checkIfExists: true)
  }
  
  if ( params.trim ) {
    (ch_reads,ch_trimgalore_txt) = trimgalore(ch_reads_input)
    ch_fastqc_trim = fastqc_trim(ch_reads)
  }
  else {
    ch_reads = ch_reads_input
    ch_fastqc = fastqc(ch_reads_input)
  }
  
  if ( params.bowtie2 ) {
    ch_align = bowtie2_align(ch_reads)
  }
  if ( params.subread ) {
    ch_align = subread_align(ch_reads)
  }
  
  ch_filter = filter(ch_align)
  
  ch_cov = coverage(ch_filter)
  
  if ( params.trim ) {
      ch_multiqc = multiqc_trim(ch_trimgalore_txt.collect(),
                                ch_fastqc_trim.collect(),
                                ch_align.collect(),
                                ch_filter.collect(),
                                ch_cov.collect())
  }
  else {
    ch_multiqc = multiqc(ch_fastqc.collect(),
                         ch_align.collect(),
                         ch_filter.collect(),
                         ch_cov.collect())
  }
  
}
