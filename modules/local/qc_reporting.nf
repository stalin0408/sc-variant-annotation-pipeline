nextflow.enable.dsl=2
//////////////////////////////////////////////////////
// PROCESS 4: QC SUMMARY REPORT
//////////////////////////////////////////////////////

process QC_SUMMARY {
    publishDir "${params.outdir}/summary",mode: 'copy', overwrite: true

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

    publishDir "${params.outdir}/combined_summary",mode: 'copy', overwrite: true


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

