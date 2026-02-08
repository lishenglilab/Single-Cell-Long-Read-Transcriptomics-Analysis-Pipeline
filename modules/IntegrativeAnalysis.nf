nextflow.enable.dsl=2

/*
 * Integrative Analysis Workflow Module
 * Combines results from multiple analyses for comprehensive interpretation
 */

process CREATE_STAGE_MAPPING_FILE {
  label 'small'
  tag "stage_mapping"
  
  publishDir "${params.outdir ?: 'Result'}/IntegrativeAnalysis/config", mode: 'copy'
  
  input:
    val stage_mapping_list
  
  output:
    path "sample_stage_mapping.tsv", emit: stage_mapping_file
  
  script:
    """
    # Create stage mapping file
    echo -e "sample\\tstage" > sample_stage_mapping.tsv
    
    # Add each mapping entry - use double quotes to allow variable expansion
    # The stage_mapping_list already contains formatted lines with sample and stage
    echo "${stage_mapping_list}" >> sample_stage_mapping.tsv
    
    echo "Stage mapping file created with \$(wc -l sample_stage_mapping.tsv) entries"
    """
}

process DIFFERENTIAL_ANALYSIS {
  label 'sclong'
  tag { comparison_name }
  
  publishDir "${params.outdir ?: 'Result'}/IntegrativeAnalysis/${comparison_name}", mode: 'copy'
  
  input:
    path seurat_tr
    path seurat_gene
    path gtf_file
    path dominant_results_dir
    path stage_mapping_file
    val ident_1
    val ident_2
    val comparison_name
  
  output:
    path "${comparison_name}/", emit: analysis_results
    path "${comparison_name}_summary.txt", emit: analysis_summary
  
  script:
    // Create output directory for this comparison
    output_dir = comparison_name
    
    // Get parameters with defaults (these should be set in main config)
    def logfc_threshold = params.integration_logfc_threshold ?: 1.5
    def pval_threshold = params.integration_pval_threshold ?: 0.05
    def min_pct = params.integration_min_pct ?: 0.01
    def n_cores = params.integration_n_cores ?: 10
    
    """
    # Create output directory
    mkdir -p ${output_dir}
    
    # Run differential analysis
    RSCRIPT_PATH="${projectDir}/bin/IntegrativeAnalysis/differential_analysis.R"
    
    Rscript \${RSCRIPT_PATH} \
      --seurat_tr_rds "${seurat_tr}" \
      --seurat_gene_rds "${seurat_gene}" \
      --gtf_file "${gtf_file}" \
      --dominant_transcripts_dir "${dominant_results_dir}" \
      --output_dir "${output_dir}" \
      --ident_1 "${ident_1}" \
      --ident_2 "${ident_2}" \
      --comparison_name "${comparison_name}" \
      --logfc_threshold ${logfc_threshold} \
      --pval_threshold ${pval_threshold} \
      --min_pct ${min_pct} \
      --n_cores ${n_cores} \
      --stage_mapping_file "${stage_mapping_file}"
    
    # Create summary file
    echo "Differential Analysis Summary" > ${comparison_name}_summary.txt
    echo "=============================" >> ${comparison_name}_summary.txt
    echo "Timestamp: \$(date)" >> ${comparison_name}_summary.txt
    echo "Comparison: ${comparison_name} (${ident_1} vs ${ident_2})" >> ${comparison_name}_summary.txt
    echo "Parameters:" >> ${comparison_name}_summary.txt
    echo "- LogFC threshold: ${logfc_threshold}" >> ${comparison_name}_summary.txt
    echo "- P-value threshold: ${pval_threshold}" >> ${comparison_name}_summary.txt
    echo "- Min percentage: ${min_pct}" >> ${comparison_name}_summary.txt
    echo "- Number of cores: ${n_cores}" >> ${comparison_name}_summary.txt
    echo "" >> ${comparison_name}_summary.txt
    echo "Output files generated in: ${output_dir}" >> ${comparison_name}_summary.txt
    echo "Generated files:" >> ${comparison_name}_summary.txt
    find ${output_dir} -type f -name "*.txt" -o -name "*.pdf" -o -name "*.png" -o -name "*.rds" | sort >> ${comparison_name}_summary.txt
    """
}

process DIFFERENTIAL_INTERPRETATION {
  label 'sclong'
  tag { comparison_name }
  
  publishDir "${params.outdir ?: 'Result'}/IntegrativeAnalysis/${comparison_name}/interpretation", mode: 'copy'
  
  // Only run when API key is provided
  when:
    params.api_key && params.api_key != '' && params.api_key != '""'
  
  input:
    path analysis_results_dir
    val comparison_name
  
  output:
    path "differential_interpretation.txt", emit: interpretation_report
    path "deepseek_differential_response.json", emit: api_response
    path "differential_interpretation_summary.txt", emit: gene_summary
  
  script:
    // Set study description
    def study_description = params.study_description 
    
    // Get API parameters
    def api_key = params.api_key ?: ''
    def model = params.deepseek_model ?: 'deepseek-reasoner'
    
    """
    # Check if API key is provided
    if [[ -z "${api_key}" ]] || [[ "${api_key}" == "" ]]; then
      echo "Differential analysis biological interpretation skipped: API key not configured" >&2
      touch differential_interpretation.txt
      echo "No API key provided, skipping DeepSeek differential interpretation" > differential_interpretation.txt
      touch deepseek_differential_response.json
      echo '{"status": "skipped", "message": "API key not provided"}' > deepseek_differential_response.json
      touch differential_genes_summary.txt
      echo "Gene summary not generated due to missing API key" > differential_genes_summary.txt
      exit 0
    fi
    
    # Check if analysis results exist
    if [[ ! -d "${analysis_results_dir}" ]]; then
      echo "ERROR: Analysis results directory not found: ${analysis_results_dir}" >&2
      exit 1
    fi
    
    # Run the DeepSeek interpretation script
    PYTHON_SCRIPT="${projectDir}/bin/deepseek/differential_interpretation.py"
    
    python3 \${PYTHON_SCRIPT} \
      --result_dir "${analysis_results_dir}" \
      --comparison_name "${comparison_name}" \
      --output_report "differential_interpretation.txt" \
      --output_api "deepseek_differential_response.json" \
      --api_key "${api_key}" \
      --study_description "${study_description}" \
      --model "${model}"
    
    # Check if report was generated
    if [[ ! -s "differential_interpretation.txt" ]]; then
      echo "WARNING: Differential interpretation report is empty or skipped" >&2
    else
      echo "Differential interpretation generated successfully"
      echo "Report preview:"
      head -20 differential_interpretation.txt
    fi
    """
}

process ORF_CLUSTER {
  label 'scmedium'
  tag { "orf_cluster" }
  
  publishDir "${params.outdir ?: 'Result'}/IntegrativeAnalysis/ORF_Cluster", mode: 'copy'
  
  input:
    path seurat_tr
    path orf_fasta
    path gtf_file
    val resolution
    val dims
  
  output:
    path "ORF_Cluster", emit: orf_cluster_results
  
  script:
    // Define output directory
    output_dir = "ORF_Cluster"
    
    // Get parameters with defaults
    def resolution_val = resolution ?: 0.4
    def dims_val = dims ?: "1:30"
    
    """
    # Create output directory
    mkdir -p ${output_dir}
    
    # Run ORF clustering analysis
    RSCRIPT_PATH="${projectDir}/bin/IntegrativeAnalysis/orf_cluster.R"
    
    # Check if R script exists
    if [ ! -f "\${RSCRIPT_PATH}" ]; then
      echo "Error: R script not found at \${RSCRIPT_PATH}" >&2
      exit 1
    fi
    
    Rscript \${RSCRIPT_PATH} \
      --seurat_rds "${seurat_tr}" \
      --orf_fasta "${orf_fasta}" \
      --gtf_file "${gtf_file}" \
      --output_dir "${output_dir}" \
      --resolution ${resolution_val} \
      --dims "${dims_val}"
    
    # Check if analysis completed successfully
    if [ ! -f "${output_dir}/sce.orf.rds" ]; then
      echo "Error: ORF clustering analysis failed - output file not generated" >&2
      exit 1
    fi
    
    # Create summary of generated files
    echo "ORF cluster analysis completed successfully!"
    echo "Generated files:"
    echo "----------------"
    find ${output_dir} -type f | sort | while read file; do
      size=\$(du -h "\$file" | cut -f1)
      echo "- \${file##*/} (\$size)"
    done
    
    # Count lines in key output files for quick verification
    echo ""
    echo "File contents summary:"
    echo "----------------------"
    if [ -f "${output_dir}/ORF_id_map.txt" ]; then
      lines=\$(wc -l < "${output_dir}/ORF_id_map.txt")
      echo "- ORF_id_map.txt: \$lines lines"
    fi
    
    if [ -f "${output_dir}/orf_aggregated_expression.csv" ]; then
      lines=\$(wc -l < "${output_dir}/orf_aggregated_expression.csv")
      echo "- orf_aggregated_expression.csv: \$lines lines"
    fi
    
    # Print cluster information from summary
    if [ -f "${output_dir}/orf_cluster_summary.txt" ]; then
      echo ""
      echo "Cluster analysis summary:"
      echo "------------------------"
      grep -A 10 "Clusters found" "${output_dir}/orf_cluster_summary.txt" || true
    fi
    """
}
workflow IntegrativeAnalysis {
  take:
    seurat_tr
    seurat_gene
    gtf_file
    dominant_results
    orf_fasta  // Add orf_fasta as input parameter
  
  main:
    // Get comparisons from config
    def comparisons = params.comparisons ?: [
      ["AD", "Normal", "AD_vs_Normal"]
    ]
    
    // Prepare stage mapping content
    def stage_mapping_content = new StringBuilder()
    params.sample_stage_mapping.each { mapping ->
      stage_mapping_content.append("${mapping[0]}\t${mapping[1]}\n")
    }
    
    // Create stage mapping file
    stage_mapping = CREATE_STAGE_MAPPING_FILE(stage_mapping_content.toString())
    
    // Run ORF clustering analysis
    orf_cluster_results = ORF_CLUSTER(
      seurat_tr,
      orf_fasta,  // Use the input parameter
      gtf_file,
      params.orf_resolution ?: 0.4,
      params.orf_dims ?: "1:30"
    )
    
    // Process each comparison
    diff_analysis_results = []
    differential_interpretation_results = []
    
    comparisons.each { comparison ->
      def ident_1 = comparison[0]
      def ident_2 = comparison[1]
      def comparison_name = comparison[2]
      
      // Run differential analysis
      diff_results = DIFFERENTIAL_ANALYSIS(
        seurat_tr,
        seurat_gene,
        gtf_file,
        dominant_results,
        stage_mapping.stage_mapping_file,
        ident_1,
        ident_2,
        comparison_name
      )
      
      diff_analysis_results.add(diff_results)
      
      // Run differential interpretation (if API key is provided)
      diff_interpretation = DIFFERENTIAL_INTERPRETATION(
        diff_results.analysis_results,
        comparison_name
      )
      
      differential_interpretation_results.add(diff_interpretation)
    }
  
  emit:
    // Stage mapping file
    stage_mapping_file = stage_mapping.stage_mapping_file
    
    // ORF clustering results
    orf_cluster = orf_cluster_results.orf_cluster_results
    
    // Differential analysis results
    diff_analysis_results = diff_analysis_results.analysis_results
}