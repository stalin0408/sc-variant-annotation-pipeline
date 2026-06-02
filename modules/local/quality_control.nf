nextflow.enable.dsl=2

/*
 * Quality Control Module for Single-cell Data
 */

//////////////////////////////////////////////////////
// PROCESS 1: Validate BAM
//////////////////////////////////////////////////////

process VALIDATE_BAM {

    label 'process_low'

    //container 'quay.io/biocontainers/samtools:1.17--h00cdaf9_0'
//    container "${params.docker_registry}/${params.docker_image}"

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path(bam), path(bai), emit: validated_bam
    tuple val(sample_id), path("${sample_id}_validation.txt"), emit: validation_report

    script:
    """
    echo "Validating BAM file: ${bam}" > ${sample_id}_validation.txt
    echo "================================" >> ${sample_id}_validation.txt

    if samtools quickcheck -v "${bam}" 2>> "${sample_id}_validation.txt"; then
        echo "PASS: BAM file is valid" >> "${sample_id}_validation.txt"
    else
        echo "FAIL: BAM file is corrupted" >> "${sample_id}_validation.txt"
    fi

    echo "" >> ${sample_id}_validation.txt
    echo "Checking index..." >> ${sample_id}_validation.txt

    if [ -f "${bai}" ]; then
        echo "PASS: BAM index exists" >> "${sample_id}_validation.txt"
    else
        echo "WARNING: BAM index missing, creating one..." >> "${sample_id}_validation.txt"
        samtools index "${bam}"
    fi

    echo "" >> ${sample_id}_validation.txt
    echo "Basic Statistics:" >> ${sample_id}_validation.txt
    samtools flagstat "${bam}" >> ${sample_id}_validation.txt
    """
}

//////////////////////////////////////////////////////
// PROCESS 2: Calculate Cell QC Metrics
//////////////////////////////////////////////////////

process CELL_QC_METRICS {

    label 'process_medium'

    publishDir "${params.outdir}/qc",mode: 'copy',overwrite: true

    //container 'quay.io/biocontainers/scanpy:1.9.1--pyhdfd78af_0'

    input:
    tuple val(sample_id), path(bam), path(barcode)

    output:
    tuple val(sample_id), path("${sample_id}_cell_metrics.tsv"), emit: cell_metrics
    tuple val(sample_id), path("${sample_id}_cell_qc.pdf"), emit: qc_plots

script:
"""

python3 << 'PYCODE'

import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sample_id = "${sample_id}"
bam_file = "${bam}"
barcode_file = "${barcode}"

print(f"Processing sample: {sample_id}")

with open(barcode_file) as f:
    valid_barcodes = set(line.strip() for line in f)

bam = pysam.AlignmentFile(bam_file, "rb")

cell_metrics = {}

for read in bam:

    if not read.has_tag('CB'):
        continue

    cell_barcode = read.get_tag('CB')

    if cell_barcode not in valid_barcodes:
        continue

    if cell_barcode not in cell_metrics:
        cell_metrics[cell_barcode] = {
            'total_reads': 0,
            'mapped_reads': 0,
            'unique_reads': 0,
            'duplicate_reads': 0
        }

    cell_metrics[cell_barcode]['total_reads'] += 1

    if not read.is_unmapped:
        cell_metrics[cell_barcode]['mapped_reads'] += 1

    if read.is_duplicate:
        cell_metrics[cell_barcode]['duplicate_reads'] += 1
    else:
        cell_metrics[cell_barcode]['unique_reads'] += 1

bam.close()

metrics_df = pd.DataFrame.from_dict(cell_metrics, orient='index')
metrics_df.index.name = 'cell_barcode'

# Safety check for empty metrics
if len(metrics_df) == 0:
    print("No valid barcodes found in BAM")
    metrics_df = pd.DataFrame(columns=[
        'total_reads',
        'mapped_reads',
        'unique_reads',
        'duplicate_reads',
        'mapping_rate',
        'duplicate_rate',
        'sample_id'
    ])
else:
    metrics_df['mapping_rate'] = metrics_df['mapped_reads'] / metrics_df['total_reads'] * 100
    metrics_df['duplicate_rate'] = metrics_df['duplicate_reads'] / metrics_df['total_reads'] * 100
    metrics_df['sample_id'] = sample_id

metrics_df.to_csv(f"{sample_id}_cell_metrics.tsv", sep="\\t")

if len(metrics_df) > 0:
    fig, axes = plt.subplots(2,2, figsize=(12,10))

    axes[0,0].hist(metrics_df['total_reads'], bins=50)
    axes[0,0].set_title("Reads per Cell")

    axes[0,1].hist(metrics_df['mapping_rate'], bins=50)
    axes[0,1].set_title("Mapping Rate")

    axes[1,0].hist(metrics_df['duplicate_rate'], bins=50)
    axes[1,0].set_title("Duplicate Rate")

    sorted_reads = np.sort(metrics_df['total_reads'])[::-1]
    cumulative = np.cumsum(sorted_reads) / np.sum(sorted_reads) * 100

    axes[1,1].plot(range(1,len(cumulative)+1), cumulative)
    axes[1,1].set_title("Cumulative Reads")

    plt.tight_layout()
    plt.savefig(f"{sample_id}_cell_qc.pdf")

print("QC metrics generated")

PYCODE
"""
}

