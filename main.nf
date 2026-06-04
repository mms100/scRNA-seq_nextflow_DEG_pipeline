nextflow.enable.dsl = 2

// --- Parameters ---
params.results_dir   = "${PWD}/results" 
params.outdir        = "${params.results_dir}/MAST_${params.cond1}vs${params.cond2}"
params.species       = "mouse"

process runMAST {
    // This goes into the raw_results/tables folder[cite: 1]
    publishDir "${params.outdir}/raw_results/tables", mode: 'copy', overwrite: true, pattern: "*.csv"
    
    input:
    val ready

    output:
    path "*.csv", emit: csv_files

    script:
    """
    Rscript ${PWD}/scripts/MAST_rcript.R \
        --object "${params.object}" \
        --cond1 ${params.cond1} \
        --cond2 ${params.cond2} \
        --cond_colname ${params.cond_colname} \
        --batch_colname ${params.batch_colname} \
        --annotation ${params.annotation} \
        --outdir .
    """
}

process filterResults {
    // This goes into the filtered_results/tables folder[cite: 8]
    publishDir "${params.outdir}/filtered_results/tables", mode: 'copy', overwrite: true
    
    input:
    path raw_csvs

    output:
    path "filtered*.csv", emit: filtered_csvs

    script:
    """
    Rscript ${PWD}/scripts/filter_script.R --input_dir . --outdir . --species "${params.species}"
    """
}

process processVolcano {
    // Published into: [raw or filtered]/plots/volcano
    publishDir "${params.outdir}/${type}/plots/volcano", mode: 'copy', overwrite: true

    input:
    tuple val(type), path(csv_files)

    output:
    path "*.pdf"

    script:
    """
    Rscript ${PWD}/scripts/volcano_plot_rscript.R \
        --input_dir_1 . \
        --outdir_3 . \
        --cond1 ${params.cond1} \
        --cond2 ${params.cond2}
    """
}

process processTop20 {
    // Published into: [raw or filtered]/plots/top_gens
    publishDir "${params.outdir}/${type}/plots/top_gens", mode: 'copy', overwrite: true

    input:
    tuple val(type), path(csv_files)

    output:
    path "*.pdf"

    script:
    """
    Rscript ${PWD}/scripts/barplot_top_20.R \
        --input_dir_1 . \
        --cond1 ${params.cond1} \
        --cond2 ${params.cond2}
    """
}

process processSummaryBarplot {
    // Published into: [raw or filtered]/plots
    publishDir "${params.outdir}/${type}/plots", mode: 'copy', overwrite: true

    input:
    tuple val(type), path(csv_files)

    output:
    path "barplot.pdf"

    script:
    """
    Rscript ${PWD}/scripts/bar_plot_rscript.R \
        --input_dir_1 . \
        --outdir .
    """
}

workflow {
    raw_ch   = runMAST(true)
    filt_ch  = filterResults(raw_ch)

    // Tagging the data streams to determine the subfolder names[cite: 1, 8]
    ch_to_plot = Channel.from("raw_results").combine(raw_ch.collect().toList())
                 .mix(Channel.from("filtered_results").combine(filt_ch.collect().toList()))

    processVolcano(ch_to_plot)
    processTop20(ch_to_plot)
    processSummaryBarplot(ch_to_plot)
}