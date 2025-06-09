// Enable DSL2 syntax
nextflow.enable.dsl = 2

// Define base directory and subdirectories
def outdir = "${PWD}/MAST_${params.cond1}vs${params.cond2}_${params.output_1}"
def tables = "${outdir}/tables"
def plots = "${outdir}/plots"
def volcano = "${outdir}/plots/volcano"
def top_gens = "${outdir}/plots/top_gens"

// Process to ensure all directories exist
process createDirectories {
    executor = 'local'
    clusterOptions = '--ntasks=1 --mem=1Gb --time=00:05:00'

    output:
    path "dummy1.txt"

    script:
    """
    mkdir -p "${outdir}"
    mkdir -p "${tables}"
    mkdir -p "${plots}"
    mkdir -p "${volcano}"
    mkdir -p "${top_gens}"
    touch dummy1.txt
    """
}

// Define a process that runs an R script
process runRscript {
    executor = 'slurm'
    clusterOptions = '--ntasks=1  --mem=45Gb --time=24:00:00'
    def rscript_path = "${PWD}/scripts/MAST_rcript_1.R"

    input:
    path dummy1

    output:
    path "dummy2.txt"
    file("*.csv")

    publishDir "${outdir}/tables", mode: 'copy'

    script:
    """
    Rscript ${rscript_path} --object "${params.object}" --cond1 ${params.cond1} --cond2 ${params.cond2} --cond_colname ${params.cond_colname} --batch_colname ${params.batch_colname} --annotation ${params.annotation} --outdir "${tables}"
    touch dummy2.txt
    """
}

process processvolcano {
    executor = 'slurm'
    clusterOptions = '--ntasks=1 --mem=10Gb --time=02:00:00'
    def rscript_3_path = "${PWD}/scripts/volcano_plot_rscript.R"

    input:
    path dummy2
    file("*.csv")

    output:
    path "dummy3.txt"
    file("*.pdf")

    publishDir volcano, mode: 'copy', overwrite: true, pattern: "*.pdf"

    script:
    """
    Rscript ${rscript_3_path} --input_dir_1 "${tables}" --outdir_3 "${volcano}" --cond1 ${params.cond1} --cond2 ${params.cond2}
    touch dummy3.txt
    """
}

process processtop_20 {
    executor = 'slurm'
    clusterOptions = '--ntasks=1 --mem=10Gb --time=02:00:00'
    def rscript_4_path = "${PWD}/scripts/barplot_top_20.R"

    input:
    path dummy3
    file("*.csv")

    output:
    path "dummy4.txt"
    file("*.pdf")

    publishDir top_gens, mode: 'copy', overwrite: true, pattern: "*.pdf"

    script:
    """
    Rscript ${rscript_4_path} --input_dir_1 "${tables}" --cond1 ${params.cond1} --cond2 ${params.cond2}
    touch dummy4.txt
    """
}

process processCSVFiles {
    executor = 'slurm'
    clusterOptions = '--ntasks=1 --mem=10Gb --time=02:00:00'
    def rscript_2_path = "${PWD}/scripts/bar_plot_rscript.R"

    input:
    path dummy4
    file("*.csv")

    output:
    file("barplot.pdf")

    publishDir plots, mode: 'copy', overwrite: true, pattern: "barplot.pdf"

    script:
    """
    Rscript ${rscript_2_path} --input_dir_1 "${tables}" --outdir "${tables}"
    """
}

workflow {
    ch1 = createDirectories()
    ch2 = runRscript(ch1.out)
    ch3 = processvolcano(ch2.out, runRscript.out)
    ch4 = processtop_20(ch3.out, runRscript.out)
    processCSVFiles(ch4.out, runRscript.out)
}
