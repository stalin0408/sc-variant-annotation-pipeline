#!/usr/bin/env nextflow

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

//////////////////////////////////////////////////////
// PROCESS 3: Filter Cells
//////////////////////////////////////////////////////

process FILTER_CELLS {

    label 'process_medium'

    input:
    tuple val(sample_id), path(cell_metrics)

    output:
    tuple val(sample_id), path("${sample_id}_filtered_cells.txt"), emit: filtered_cells
    tuple val(sample_id), path("${sample_id}_filtered_metrics.tsv"), emit: filtered_metrics

    script:
"""
python3 << 'PYCODE'

import pandas as pd

df = pd.read_csv("${cell_metrics}", sep="\\t", index_col=0)

required_cols = ["total_reads", "mapping_rate", "duplicate_rate"]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

min_reads = ${params.min_reads_per_cell ?: 1000}
min_mapping = ${params.min_mapping_rate ?: 50}
max_dup = ${params.max_duplicate_rate ?: 50}

filtered = df[
    (df.total_reads >= min_reads) &
    (df.mapping_rate >= min_mapping) &
    (df.duplicate_rate <= max_dup)
]

print(f"Input cells: {len(df)}")
print(f"Filtered cells: {len(filtered)}")

filtered.index.to_series().to_csv("${sample_id}_filtered_cells.txt", header=False)
filtered.to_csv("${sample_id}_filtered_metrics.tsv", sep="\\t")

PYCODE

"""
}

//////////////////////////////////////////////////////
// PROCESS 4: QC SUMMARY REPORT
//////////////////////////////////////////////////////

process QC_SUMMARY {

    label 'process_low'

    input:
    tuple val(sample_id), path(validation), path(cell_metrics), path(filtered_cells)

    output:
    tuple val(sample_id), path("${sample_id}_qc_summary.html"), emit: qc_summary

    script:
    """
    python3 << 'PYCODE'
import pandas as pd
from datetime import datetime
import html

metrics = pd.read_csv("${cell_metrics}", sep="\t", index_col=0)

# ✅ FIX validation file handling
try:
    with open("${validation}", "r", encoding="utf-8", errors="replace") as f:
        validation_text = html.escape(f.read())
except Exception as e:
    validation_text = f"Error reading validation file: {e}"

total_cells = len(metrics)

# ✅ FIX filtered_cells handling
try:
    with open("${filtered_cells}") as f:
        passed = sum(1 for _ in f)
except Exception as e:
    print(f"Error reading filtered cells file: {e}")
    passed = 0

# ✅ ALWAYS define html_content (outside try)
html_content = f\"\"\"
<html>
<h1>QC Summary - ${sample_id}</h1>
<p>Generated: {datetime.now()}</p>

<h2>Cells</h2>
Total: {total_cells}<br>
Passed: {passed}

<h2>BAM Validation</h2>
<pre>{validation_text}</pre>

</html>
\"\"\"

with open("${sample_id}_qc_summary.html","w") as f:
    f.write(html_content)
PYCODE
"""
}

//////////////////////////////////////////////////////
// PROCESS 5: Combine Metrics Across Samples
//////////////////////////////////////////////////////
process COMBINE_METRICS {

    label 'process_low'

    input:
    path metrics_files

    output:
    path "all_samples_filtered_metrics.tsv"

    script:
    """
    echo "==== FILES IN WORKDIR ===="
    ls -lh
    python3 - << 'PYCODE'

import pandas as pd
import glob
import os

output_file = "all_samples_filtered_metrics.tsv"
files = sorted(
    f for f in glob.glob("*_filtered_metrics.tsv")
    if f != output_file
)

for f in files:
    print(f, "exists?", os.path.exists(f))

print("DEBUG files:", files)

if len(files) == 0:
    raise ValueError("No metric files provided")

dfs = [pd.read_csv(f, sep="\\t", index_col=0) for f in files]

combined = pd.concat(dfs, axis=0)

combined.to_csv(output_file, sep="\\t")

print(f"Combined {len(files)} files")

PYCODE
"""
}


//////////////////////////////////////////////////////
// WORKFLOW
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
