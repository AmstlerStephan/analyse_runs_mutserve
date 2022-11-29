process SUMMARIZE_NGS_UMI {
  input:
    tuple val( run ), path( ngs_data )
    path summarize_results_R
    val data_type
    val output_type
  output:
    path 
  script:
  """
  
  """
}