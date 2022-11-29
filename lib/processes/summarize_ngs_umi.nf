process SUMMARIZE_NGS_UMI {
    publishDir "${params.output}/NGS/", mode: 'symlink'
  input:
    tuple val( run ), path( ngs_data )
    path summarize_results_R
    val data_type
    val output_type
  output:
    path ("umi*")
  script:
  """
  
  """
}