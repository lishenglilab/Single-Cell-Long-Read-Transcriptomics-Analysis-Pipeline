// main.nf
nextflow.enable.dsl=2

include { InputAndPreprocess } from './modules/InputAndPreprocess.nf'
include { AlignmentAndTranscriptReconstruction } from './modules/AlignmentAndTranscriptReconstruction.nf'
include { PreliminarySingleCellAnalysis } from './modules/PreliminarySingleCellAnalysis.nf'
include { AlternativeSplicing } from './modules/AlternativeSplicing.nf'
include { IntegrativeAnalysis } from './modules/IntegrativeAnalysis.nf'

workflow {
  // Module (1): InputAndPreprocess
  pre = InputAndPreprocess(params.metadata, params.whitelist)

  // Module (2): AlignmentAndTranscriptReconstruction
  align = AlignmentAndTranscriptReconstruction(
    pre.filtered_fastq,
    params.reference_fa,
    params.genedb_gtf
  )

  // Module (3): PreliminarySingleCellAnalysis
  sc = PreliminarySingleCellAnalysis(align.isoquant_out)

  // Module (4): AlternativeSplicing
  as = AlternativeSplicing(
    sc.seurat_tr,
    sc.seurat_gene,
    align.models_gtf,         
    align.models_fa_1line,    
    align.isoquant_out        
  )
  
  // Module (5): IntegrativeAnalysis 
  integrate = IntegrativeAnalysis(
    sc.seurat_tr,           // seurat_tr
    sc.seurat_gene,         // seurat_gene
    align.models_gtf,       // gtf_file
    as.dominant_results,    // dominant_results
    align.orf_fa            // orf_fasta
  )

  integration_results = integrate.diff_analysis_results
  orf_cluster_results = integrate.orf_cluster


}