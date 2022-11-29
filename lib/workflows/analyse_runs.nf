#!/usr/bin/env nextflow

nextflow.enable.dsl=2

requiredParams = [
    'run_folder', 'nanostat_folder', 'nanostat_tsv_pattern', 'mutserve_summary_pattern', 'ngs_data', 'output'
]

for (param in requiredParams) {
    if (params[param] == null) {
      exit 1, "Parameter ${param} is required."
    }
}

// DEFINE PATHS # these are strings which are used to define input Channels,
// but they are specified here as they may be referenced in LOGGING
merge_result_files =                            file("../bin/merge_result_files.R", checkIfExists: true)
summarize_UMI_run_mutserve_ngs_data =           file("../bin/summarize_UMI_run_mutserve_ngs_data.R", checkIfExists: true)
summarize_results =                             file("../bin/summarize_results.R", checkIfExists: true)
num_of_read_vs_bam_files =                      file("../bin/num_of_read_vs_bam_files.R", checkIfExists: true)
summarize_results_plasmid_quality_measures =    file("../bin/summarize_results_Plasmid_quality_measures.R", checkIfExists: true)
summarize_results_plasmid_plots =               file("../bin/summarize_results_Plasmid_plots.R", checkIfExists: true)
summarize_results_NGS_quality_measures =        file("../bin/summarize_results_NGS_quality_measures.R", checkIfExists: true)
summarize_results_NGS_plots =                   file("../bin/summarize_results_NGS_plots.R", checkIfExists: true)

ngs_data = file("${params.ngs_data}", checkIfExists: true)
expected_mutations = file("${params.expected_mutations}", checkIfExists: true)

// Set directory names for publishing data
raw = "raw"
filtered = "filtered"
plots = "plots"
qm = "quality measures"


// STAGE CHANNELS
if (params.all_runs) {
    nanostat_summary_files = Channel.fromPath("${params.nanostat_folder}/run*/*${nanostat_tsv}", type: 'file')
    mutserve_summary_files = Channel.fromPath("${params.run_folder}/run*/**${params.mutserve_summary}", type: 'file')
}else{
    nanostat_summary_files = Channel.fromPath("${params.nanostat_folder}/*${nanostat_tsv}", type: 'file')
    mutserve_summary_files = Channel.fromPath("${params.run_folder}/**${params.mutserve_summary}", type: 'file')
}

nanostat_summary_files
.map{ 
    nanostat_summary_path ->
    run = (nanostat_summary_path =~ /run\d*_*V*\d*/)[0]
    tuple( run, nanostat_summary_path )
}
.set{ nanostat_summaries }

mutserve_summary_files
.map{
    mutserve_summary_path ->
    run = (nanostat_summary_path =~ /run\d*_*V*\d*/)[0]
    tuple( run, mutserve_summary_path )
}
.join(nanostat_summaries)
.view()
.set{ run_summaries }

include {SUMMARIZE_RUN} from '../processes/summarize_run.nf'
include {MERGE_RESULT_FILES} from '../processes/merge_result_files.nf'
include {SUMMARIZE_RESULTS} from '../processes/summarize_results.nf'
include {SUMMARIZE_NGS_UMI as SUMMARIZE_NGS_UMI_PLOTS; 
        SUMMARIZE_NGS_UMI as SUMMARIZE_NGS_UMI_QM;
        SUMMARIZE_NGS_UMI as SUMMARIZE_NGS_UMI_FILTERED_PLOTS; 
        SUMMARIZE_NGS_UMI as SUMMARIZE_NGS_UMI_FILTERED_QM;
        } from '../processes/summarize_ngs_umi.nf'
include {SUMMARIZE_PLASMID_UMI as SUMMARIZE_PLASMID_UMI_PLOTS; 
        SUMMARIZE_PLASMID_UMI as SUMMARIZE_PLASMID_UMI_QM;
        SUMMARIZE_PLASMID_UMI as SUMMARIZE_PLASMID_UMI_FILTERED_PLOTS; 
        SUMMARIZE_PLASMID_UMI as SUMMARIZE_PLASMID_UMI_FILTERED_QM;
        } from '../processes/summarize_plasmid_umi.nf'

workflow ANALYSE_RUN {

    SUMMARIZE_RUN( run_summaries, summarize_UMI_run_mutserve_ngs_data )

    if(params.merge_result_files){
        SUMMARIZE_RUN.out.umi_all_tsv
        .map{ run, umi_all_tsv -> umi_all_tsv}
        .collect()
        .set{ result_files }

        MERGE_RESULT_FILES( result_files, merge_result_files)
    }

    SUMMARIZE_NGS_UMI_PLOTS( SUMMARIZE_RUN.out.ngs_raw, summarize_results_NGS_plots, plots, raw )
    SUMMARIZE_NGS_UMI_QM( SUMMARIZE_RUN.out.ngs_raw, summarize_results_NGS_quality_measures, qm, raw )
    SUMMARIZE_NGS_UMI_FILTERED_PLOTS( SUMMARIZE_RUN.out.ngs_filtered, summarize_results_NGS_plots, plots, filtered )
    SUMMARIZE_NGS_UMI_FILTERED_QM( SUMMARIZE_RUN.out.ngs_filtered, summarize_results_NGS_quality_measures, qm, filtered )
    SUMMARIZE_PLASMID_UMI_PLOTS( SUMMARIZE_RUN.out.plasmids_raw, summarize_results_plasmid_plots, plots, raw ) 
    SUMMARIZE_PLASMID_UMI_QM( SUMMARIZE_RUN.out.plasmids_raw, summarize_results_plasmid_quality_measures, qm, raw )
    SUMMARIZE_PLASMID_UMI_FILTERED_PLOTS( SUMMARIZE_RUN.out.plasmids_raw, summarize_results_plasmid_plots, plots, filtered )
    SUMMARIZE_PLASMID_UMI_FILTERED_QM( SUMMARIZE_RUN.out.plasmids_raw, summarize_results_plasmid_quality_measures, qm, filtered )




}