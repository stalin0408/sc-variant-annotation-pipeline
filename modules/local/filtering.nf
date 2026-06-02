nextflow.enable.dsl=2
//////////////////////////////////////////////////////
// PROCESS 3: Filter Cells
//////////////////////////////////////////////////////

process FILTER_CELLS {
    publishDir "${params.outdir}/filtered_cells",mode: 'copy'

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

