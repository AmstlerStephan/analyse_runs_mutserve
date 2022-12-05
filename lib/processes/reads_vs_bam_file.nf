process READS_VS_BAM_FILE {
  publishDir "${params.output}/reads_vs_bam_file/", mode: 'copy'
  input:
    path umi_mutserve_merged
    path num_of_reads_vs_bam_file_R
  output:
    path ("umi*")
  script:
  """
    Rscript ${num_of_reads_vs_bam_file_R} --umi_merged ${umi_mutserve_merged} 
  """
}