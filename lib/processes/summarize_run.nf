process SUMMARIZE_RUN {
    publishDir "${params.output}/umi_summary_files_per_run/${run}/", mode: 'copy'
  input:
    tuple val( run ), path( mutserve_summary ), path( nanostat_summary )
    path corresponding_position
    path ngs_data
    path summarize_UMI_run_mutserve_ngs_data_R
  output:
    tuple val( "${run}"), path( "UMI_sequencing_mutserve_all*"), emit: umi_all
    tuple val( "${run}"), path( "UMI_sequencing_mutserve_plasmids.tsv"), emit: plasmids_raw
    tuple val( "${run}"), path( "UMI_sequencing_mutserve_plasmids_filtered.tsv"), emit: plasmids_filtered
    tuple val( "${run}"), path( "NGS_UMI_samples.tsv"), emit: ngs_raw
    tuple val( "${run}"), path( "NGS_UMI_samples_filtered.tsv"), emit: ngs_filtered
  script:
  """
    Rscript ${summarize_UMI_run_mutserve_ngs_data_R} \
    --run ${run} \
    --nanostat_summary ${nanostat_summary} \
    --mutserve_summary ${mutserve_summary} \
    --ngs_data ${ngs_data} \
    --umi_cutoff_R9 ${params.umi_cutoff_R9} \
    --umi_cutoff_V14 ${params.umi_cutoff_V14} \
    --corresponding_position ${params.corresponding_position}
  """
}