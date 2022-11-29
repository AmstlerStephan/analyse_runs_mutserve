process  MERGE_RESULT_FILES {
    publishDir "${params.output}/umi_all/", mode: 'copy'
  input:
    path umi_sequencing_mutserve
    path merge_result_files_R
  output:
    path ("UMI_sequencing_mutserve_merged.tsv")
  script:
  """
    Rscript ${merge_result_files_R} --parsed_files ${umi_sequencing_mutserve}
  """
}