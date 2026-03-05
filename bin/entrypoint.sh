#!/bin/bash
set -e

echo "======================================"
echo "Single-cell Variant Annotation Pipeline"
echo "======================================"
echo ""

# Activate virtual environment if needed
# source /opt/venv/bin/activate

# Execute the command passed to docker run
exec "$@"

