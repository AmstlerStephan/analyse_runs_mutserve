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
merge_result_files =                            file("${projectDir}/bin/merge_result_files.R", checkIfExists: true)
summarize_UMI_run_mutserve_ngs_data =           file("${projectDir}/bin/summarize_UMI_run_mutserve_ngs_data.R", checkIfExists: true)
num_of_read_vs_bam_files =                      file("${projectDir}/bin/num_of_read_vs_bam_files.R", checkIfExists: true)
summarize_results_plasmid_quality_measures =    file("${projectDir}/bin/summarize_results_Plasmid_quality_measures.R", checkIfExists: true)
summarize_results_plasmid_plots =               file("${projectDir}/bin/summarize_results_Plasmid_plots.R", checkIfExists: true)
summarize_results_NGS_quality_measures =        file("${projectDir}/bin/summarize_results_NGS_quality_measures.R", checkIfExists: true)
summarize_results_NGS_plots =                   file("${projectDir}/bin/summarize_results_NGS_plots.R", checkIfExists: true)

ngs_data =                                      file("${params.ngs_data}", checkIfExists: true)
mutation_classification =                       file("${params.expected_mutations}", checkIfExists: true)
corresponding_positions =                       file("${params.corresponding_positions}", checkIfExists: true)

// Set directory names for publishing data
raw = "raw"
filtered = "filtered"
plots = "plots"
qm = "quality measures"

print params.all_runs
print params.mutserve_dir

// STAGE CHANNELS
if (params.all_runs) {
    Channel.fromPath("${params.run_folder}/run*/${params.ont_pl_folder}/${params.mutserve_summary_pattern}", type: 'file', maxDepth: 2)
    .view()
    .set{ mutserve_summary_files}
    
    Channel.fromPath("${params.nanostat_folder}/run*/*${params.nanostat_tsv_pattern}", type: 'file')
    .view()
    .set{ nanostat_summary_files }
}else if (params.mutserve_dir) {
    Channel.fromPath("${params.run_folder}/**${params.mutserve_summary_pattern}", type: 'file')
    .set{ mutserve_summary_files}
    
    Channel.fromPath("${params.nanostat_folder}/run*/*${params.nanostat_tsv_pattern}", type: 'file')
    .set{ nanostat_summary_files }
}else{
    Channel.fromPath("${params.run_folder}/${params.ont_pl_folder}/**${params.mutserve_summary_pattern}", type: 'file')
    .view()
    .set{ mutserve_summary_files}
    
    Channel.fromPath("${params.nanostat_folder}/**${params.nanostat_tsv_pattern}", type: 'file')
    .view()
    .set{ nanostat_summary_files }
}

mutserve_summary_files
.map{
    mutserve_summary_path ->
    run = (mutserve_summary_path =~ /run\d*_*V*\d*/)[0]
    tuple( run, mutserve_summary_path )
}
.set{ mutserve_summaries }

nanostat_summary_files
.map{ 
    nanostat_summary_path ->
    run = (nanostat_summary_path =~ /run\d*_*V*\d*/)[0]
    tuple( run, nanostat_summary_path )
}
.set{ nanostat_summaries }

mutserve_summaries
.join( nanostat_summaries )
.set{ run_summaries }



include {SUMMARIZE_RUN} from '../processes/summarize_run.nf'
include {MERGE_RESULT_FILES} from '../processes/merge_result_files.nf'
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
include {READS_VS_BAM_FILE} from '../processes/reads_vs_bam_file.nf'

workflow ANALYSE_RUN {

    SUMMARIZE_RUN( run_summaries, ngs_data, corresponding_positions, summarize_UMI_run_mutserve_ngs_data )

    if(params.merge_umi_result_file){
        SUMMARIZE_RUN.out.umi_all
        .map{ run, umi_all_tsv -> umi_all_tsv}
        .collect()
        .set{ result_files }

        MERGE_RESULT_FILES( result_files, merge_result_files)
        READS_VS_BAM_FILE( MERGE_RESULT_FILES.out.umi_mutserve_merged, num_of_read_vs_bam_files )
    }

    SUMMARIZE_NGS_UMI_PLOTS( SUMMARIZE_RUN.out.ngs_raw, summarize_results_NGS_plots, plots, raw )
    SUMMARIZE_NGS_UMI_QM( SUMMARIZE_RUN.out.ngs_raw, summarize_results_NGS_quality_measures, qm, raw )
    SUMMARIZE_NGS_UMI_FILTERED_PLOTS( SUMMARIZE_RUN.out.ngs_filtered, summarize_results_NGS_plots, plots, filtered )
    SUMMARIZE_NGS_UMI_FILTERED_QM( SUMMARIZE_RUN.out.ngs_filtered, summarize_results_NGS_quality_measures, qm, filtered )
    SUMMARIZE_PLASMID_UMI_PLOTS( SUMMARIZE_RUN.out.plasmids_raw, mutation_classification, summarize_results_plasmid_plots, plots, raw ) 
    SUMMARIZE_PLASMID_UMI_QM( SUMMARIZE_RUN.out.plasmids_raw, mutation_classification, summarize_results_plasmid_quality_measures, qm, raw )
    SUMMARIZE_PLASMID_UMI_FILTERED_PLOTS( SUMMARIZE_RUN.out.plasmids_filtered, mutation_classification, summarize_results_plasmid_plots, plots, filtered )
    SUMMARIZE_PLASMID_UMI_FILTERED_QM( SUMMARIZE_RUN.out.plasmids_filtered, mutation_classification, summarize_results_plasmid_quality_measures, qm, filtered )

}
