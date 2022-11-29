process SUMMARIZE_PLASMID_UMI {
    publishDir "${params.output}/plasmid/${data_type}/${output_type}", mode: 'symlink'
  input:
    tuple val( run ), path( plasmid_data )
    path mutation_calssification
    path summarize_results_R
    val data_type
    val output_type
  output:
    path ("umi*")
  script:
  """
    Rscript ${summarize_results_R} --umi_plasmid_samples ${ngs_data} --mutation_classification ${mutation_calssification}
  """
}