process SUMMARIZE_NGS_UMI {
    publishDir "${params.output}/NGS/${data_type}/${output_type}/${run}", mode: 'copy'
  input:
    tuple val( run ), path( ngs_data )
    path summarize_results_R
    val data_type
    val output_type
  output:
    path ("umi*")
  script:
  """
    Rscript ${summarize_results_R} --ngs_umi_samples ${ngs_data}
  """
}