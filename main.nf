nextflow.enable.dsl = 2

// Use params for all directory paths
params.outdir   = "${PWD}/MAST_${params.cond1}vs${params.cond2}"
params.tables   = "${params.outdir}/tables"
params.plots    = "${params.outdir}/plots"
params.volcano  = "${params.outdir}/plots/volcano"
params.top_gens = "${params.outdir}/plots/top_gens"

process createDirectories {
    output:
    tuple path("dummy1.txt"), val(true)

    script:
    """
    mkdir -p "${params.outdir}"
    mkdir -p "${params.tables}"
    mkdir -p "${params.plots}"
    mkdir -p "${params.volcano}"
    mkdir -p "${params.top_gens}"
    touch dummy1.txt
    """
}

process runRscript {
    publishDir "${params.tables}", mode: 'copy',overwrite: true, pattern: "*.csv"
    input:
    tuple path(dummy1), val(x)

    output:
    tuple path("dummy2.txt"), file("*.csv")

    script:
    """
    Rscript ${PWD}/scripts/MAST_rcript.R --object "${params.object}" --cond1 ${params.cond1} --cond2 ${params.cond2} --cond_colname ${params.cond_colname} --batch_colname ${params.batch_colname} --annotation ${params.annotation} --outdir "${params.tables}"
    touch dummy2.txt
    """
}

process processvolcano {
    
    publishDir "${params.volcano}", mode: 'copy',overwrite: true, pattern: "*.pdf"

    input:
    tuple path(dummy2), file(csv_files)

    output:
    tuple path("dummy3.txt"), file("*.pdf")

    script:
    """
    Rscript ${PWD}/scripts/volcano_plot_rscript.R --input_dir_1 "${params.tables}" --outdir_3 "${params.volcano}" --cond1 ${params.cond1} --cond2 ${params.cond2}
    touch dummy3.txt
    """
}

process processtop_20 {
    
    publishDir "${params.top_gens}", mode: 'copy',overwrite: true, pattern: "*.pdf"

    input:
    tuple path(dummy3), file(csv_files)

    output:
    tuple path("dummy4.txt"), file("*.pdf")

    script:
    """
    Rscript ${PWD}/scripts/barplot_top_20.R --input_dir_1 "${params.tables}" --cond1 ${params.cond1} --cond2 ${params.cond2}
    touch dummy4.txt
    """
}

process processCSVFiles {
    
    publishDir "${params.plots}", mode: 'copy',overwrite: true, pattern: "barplot.pdf"

    input:
    tuple path(dummy4), file(csv_files)

    output:
    file("barplot.pdf")

    script:
    """
    Rscript ${PWD}/scripts/bar_plot_rscript.R --input_dir_1 "${params.tables}" --outdir "${params.tables}"
    """
}

workflow {
    ch1 = createDirectories()
    ch2 = runRscript(ch1)
    ch3 = processvolcano(ch2)
    ch4 = processtop_20(ch3)
    processCSVFiles(ch4)
}
