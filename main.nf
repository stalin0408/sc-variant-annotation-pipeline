#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


/*--------------------------------------
Single-cell Variant Annotation Pipeline
--------------------------------------*/

//------------------------------------
//Sub-workflow
//------------------------------------
include { QUALITY_CONTROL } from './subworkflows/preprocessing'
// Pipeline header



// Check required parameters
def validateParameters() {
    if (!params.input) {
    error "Input parameter is required. Please provide --input samplesheet.csv"
    }

    if (!params.outdir) {
        params.outdir = "results"
    }
}

workflow {
    log.info """
    ======================================
    Single-cell Variant Annotation Pipeline
    ======================================
    Version: 1.0.0
    """
    log.info "Parameters:"
    log.info " Input: ${params.input}"
    log.info " Output: ${params.outdir}"
    log.info " Genome: ${params.genome}"
    /*
    -------------------------
    Load sample sheet
    -------------------------
    */
    validateParameters()

    samples = Channel.fromPath(params.input)
        | splitCsv(header: true, sep: ',')
        | map { row ->

            def sample_id = row.sample_id
            def bam = file(row.bam, checkIfExists: true)
            def barcodes = file(row.barcodes, checkIfExists: true)
            def bai = row.bai ? file(row.bai, checkIfExists: true) : null

            log.info "Loaded sample: ${sample_id}"

            tuple(sample_id, bam, barcodes, bai)
        }

    /*
    -------------------------
    Run Quality Control
    -------------------------
    */

    qc_results = QUALITY_CONTROL(samples)

    /*
    -------------------------
    Save QC Outputs
    -------------------------
    */

    qc_results.qc_summaries
        | map { sample_id, html ->
            def target = file("${params.outdir}/qc/${sample_id}_qc_summary.html")
            target.parent.mkdirs()
            file(html).copyTo(target)
            tuple(sample_id, html)
        }
        | set { saved_summaries }

    qc_results.qc_plots
        | map { sample_id, pdf ->
            def target = file("${params.outdir}/qc/${sample_id}_cell_qc.pdf")
            target.parent.mkdirs()
            file(pdf).copyTo(target)
            tuple(sample_id, pdf)
        }
        | set { saved_plots }

    qc_results.filtered_cells
        | map { sample_id, txt ->
            def target = file("${params.outdir}/filtered_cells/${sample_id}_filtered_cells.txt")
            target.parent.mkdirs()
            file(txt).copyTo(target)
            tuple(sample_id, txt)
        }
        | set { saved_cells }

    qc_results.filtered_metrics
        | map { sample_id, metrics ->
            def target = file("${params.outdir}/qc/${sample_id}_filtered_metrics.tsv")
            target.parent.mkdirs()
            file(metrics).copyTo(target)
            tuple(sample_id, metrics)
        }
        | set { saved_metrics }

    /*
    -------------------------
    Prepare samples for next stage
    -------------------------
    */

    passing_samples = qc_results.validated_bams
        .join(qc_results.filtered_cells)
        .map { sample_id, bam, bai, filtered_cells ->
            tuple(sample_id, bam, bai, filtered_cells)
        }

    // Workflow completion
workflow.onComplete {

    def outdir = workflow.params?.outdir ?: "results"

    log.info "======================================"
    log.info "Pipeline Completed"
    log.info "======================================"

    if (workflow.success) {
        log.info "Status: SUCCESS ✅"
    } else {
        log.error "Status: FAILED ❌"
        }
    log.info "Results:"
    log.info "  - QC summaries: ${outdir}/qc/"
    log.info "  - Filtered cells: ${outdir}/filtered_cells/"
    log.info "  - Combined metrics: ${outdir}/all_samples_filtered_metrics.tsv"
    log.info "  - Passing samples list: ${outdir}/passing_samples.csv"

    log.info ""
    log.info "Duration: ${workflow.duration}"
    log.info "Exit status: ${workflow.exitStatus}"
    }

}
