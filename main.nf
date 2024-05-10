// Enable DSL2 syntax
nextflow.enable.dsl = 2


// Define base directory and subdirectories

def outdir = "${PWD}/MAST_${params.cond1}vs${params.cond2}"
def tables = "${outdir}/tables"
def plots = "${outdir}/plots"
def volcano = "${outdir}/plots/volcano"
def top_gens = "${outdir}/plots/top_gens"


// Process to ensure all directories exist
process createDirectories {
    executor = 'local'  // Can run on a local executor since it's a simple task
    clusterOptions = '--ntasks=1 --mem=1Gb --time=00:05:00'  // Minimal resources

    script:
    """
    # Create directories with the recursive option
    mkdir -p "${outdir}"
    mkdir -p "${tables}"
    mkdir -p "${plots}"
    mkdir -p "${volcano}"
    mkdir -p "${top_gens}"
    """
}

// Define a process that runs an R script
process runRscript {
    executor = 'slurm'

    clusterOptions = '--ntasks=1  --mem=45Gb --time=24:00:00'

    def rscript_path = "${PWD}/scripts/MAST_rcript.R"  // Correct path to the R script

    // Define the publish directory
    publishDir "${outdir}/tables", mode: 'copy'

    // Define input
    input:
    val params.object  // Ensure this is correct (single value input)

    // Define output
    output:

    file("*.csv")  // Specify the expected output

    // The script block for execution
    script:
    """
    Rscript ${rscript_path} --object "${params.object}" --cond1 ${params.cond1} --cond2 ${params.cond2} --cond_colname ${params.cond_colname}  --annotation ${params.annotation}  --outdir "${tables}"    
    """
}


// New process that will use the CSV files from the previous process
process processCSVFiles {
    executor = 'slurm'
    clusterOptions = '--ntasks=1 --mem=10Gb --time=02:00:00'

    def rscript_2_path = "${PWD}/scripts/bar_plot_rscript.R"  // Correct path to the R script
   
    publishDir plots, mode: 'copy', overwrite: true, pattern: "barplot.pdf"

    input:
    file("*.csv")   // Receiving the CSV files

    output:
    file("barplot.pdf")  // Resulting file after processing

    script:
    """
    Rscript ${rscript_2_path} --input_dir_1 "${tables}" --outdir "${tables}"

    """
}

// Another process that will use the CSV files from the previous process
process processvolcano {
    executor = 'slurm'
    clusterOptions = '--ntasks=1 --mem=10Gb --time=02:00:00'

    def rscript_3_path = "${PWD}/scripts/volcano_plot_rscript.R"  // Correct path to the R script
   
    publishDir volcano, mode: 'copy', overwrite: true, pattern: "*.pdf"

    input:
    file("*.csv")   // Receiving the CSV files

    output:
    file("*.pdf")  // Resulting file after processing

    script:
    """
    Rscript ${rscript_3_path} --input_dir_1 "${tables}" --outdir_3 "${volcano}"  --cond1 ${params.cond1} --cond2 ${params.cond2}

    """
}

// Another process that will use the CSV files from the previous process
process processtop_20 {
    executor = 'slurm'
    clusterOptions = '--ntasks=1 --mem=10Gb --time=02:00:00'

    def rscript_4_path = "${PWD}/scripts/barplot_top_20.R"  // Correct path to the R script
   
    
    publishDir top_gens, mode: 'copy', overwrite: true, pattern: "*.pdf"

    input:
    file("*.csv")   // Receiving the CSV files

    output:
    file("*.pdf")  // Resulting file after processing

    script:
    """
    Rscript ${rscript_4_path} --input_dir_1 "${tables}"   --cond1 ${params.cond1} --cond2 ${params.cond2}

    """
}



workflow {
    // Step 1: Create required directories
    createDirectories()

    // Step 2: Run R script and create a channel for CSV files
    runRscript(params.object)


    // Step 4: Process the CSV files
    processCSVFiles(runRscript.out)

    // Step 5: Process additional data (like a volcano plot)
    processvolcano(runRscript.out)

    // Step 6: Process the top 20 data
    processtop_20(runRscript.out)

}
