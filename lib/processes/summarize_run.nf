process SUMMARIZE_RUN {
    publishDir "${params.output}/UMI_vs_NGS_mutserve/${run}/*tsv", mode: 'copy'
  input:
    tuple path( run ), path( mutserve_summary ), path( nanostat_summary )
    val ngs_data
    val summarize_UMI_run_mutserve_ngs_data
  output:
    tuple path( "${run}"), path( "UMI_sequencing_mutserve.tsv"), emit: umi_all
    tuple path( "${run}"), path( "UMI_sequencing_mutserve_plasmids.tsv"), emit: plasmids
    tuple path( "${run}"), path( "UMI_sequencing_mutserve_plasmids_filtered.tsv"), emit: plasmids_filtered
    tuple path( "${run}"), path( "NGS_UMI_samples.tsv"), emit: ngs
    tuple path( "${run}"), path( "NGS_UMI_samples_filtered.tsv"), emit: ngs_filtered
  script:
  """
    Rscript ${summarize_UMI_run_mutserve_ngs_data} \
    --run ${run} \
    --nanostat_summary ${nanostat_summary} \
    --mutserve_summary ${mutserve_summary} \
    --ngs_data ${ngs_data}
  """
}