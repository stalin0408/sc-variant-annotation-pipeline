nextflow.enable.dsl=2

include { VALIDATE_BAM; CELL_QC_METRICS } from '../modules/local/quality_control.nf'
include { FILTER_CELLS } from '../modules/local/filtering.nf'
include { QC_SUMMARY; COMBINE_METRICS } from '../modules/local/qc_reporting.nf'


//////////////////////////////////////////////////////
// PREPROCESING WORKFLOW
//////////////////////////////////////////////////////

workflow QUALITY_CONTROL {

    take:
    samples  // [sample_id, bam, barcodes, bai]

    main:

    validated = VALIDATE_BAM(
        samples.map { sample_id, bam, barcodes, bai ->
            tuple(sample_id, bam, bai)
        }
    )

    metrics = CELL_QC_METRICS(
        samples.map { sample_id, bam, barcodes, bai ->
            tuple(sample_id, bam, barcodes)
        }
    )

    // ✅ Use emit-based style (recommended)
    filtered = FILTER_CELLS(
        metrics.cell_metrics.map { sample_id, file ->
            tuple(sample_id, file)
        }
    )

    // ✅ Debug
    filtered.filtered_cells.view { it -> "DEBUG FILTER: $it" }
    filtered.filtered_metrics.view { it -> "DEBUG FILTER: $it" }
    
    // ✅ Correct joins
    qc_input = validated.validation_report
        .join(metrics.cell_metrics, by: 0)
        .join(filtered.filtered_cells, by: 0)

    qc_input.view()

    summaries = QC_SUMMARY(qc_input)

    combine_input = filtered.filtered_metrics

    // ✅ Use filtered metrics (not raw)
    COMBINE_METRICS(
    combine_input
        .map { sample_id, file -> file }
    )

    emit:
    validated_bams   = validated.validated_bam
    cell_metrics     = metrics.cell_metrics
    filtered_cells   = filtered.filtered_cells
    filtered_metrics = filtered.filtered_metrics
    qc_summaries     = summaries.qc_summary
    qc_plots         = metrics.qc_plots
}
