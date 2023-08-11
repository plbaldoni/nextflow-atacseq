params.reads = "$launchDir/*{R1,R2}_001.fastq.gz"
params.outdir = "$launchDir"
params.erate = 0.2
params.mapq = 30
params.gsize = 2654621783
params.maxmultimap = 4
params.maxfraglen = 2000
params.subreadIindex = "/stornext/Home/data/allstaff/b/baldoni.p/lab_smyth/baldoni.p/software/RsubreadIndex/GRCm39/GRCm39.genome_index"
params.bowtie2Index = "/stornext/Home/data/allstaff/b/baldoni.p/lab_smyth/baldoni.p/software/bowtie2Index/GRCm39/GRCm39.genome_index"

include { fastqc } from './modules/fastqc'
include { fastqc as fastqc_trim} from './modules/fastqc'
include { trimgalore } from './modules/trimgalore'
include { bowtie2_align } from './modules/bowtie2_align'
include { filter } from './modules/filter'
include { coverage } from './modules/coverage'
include { multiqc } from './modules/multiqc'

workflow {
  
  reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
  
  ch_fastqc = fastqc(reads)
  
  (ch_trimgalore_reads,ch_trimgalore_txt) = trimgalore(reads)

  ch_fastqc_trim = fastqc_trim(ch_trimgalore_reads)
  
  ch_subread_align = bowtie2_align(ch_trimgalore_reads)
  
  ch_filter = filter(ch_subread_align)
  
  ch_cov = coverage(ch_filter)
  
  ch_multiqc = multiqc(ch_fastqc.collect(),
                       ch_trimgalore_txt.collect(),
                       ch_fastqc_trim.collect(),
                       ch_subread_align.collect(),
                       ch_filter.collect(),
                       ch_cov.collect())
}
