# Single-cell Variant Annotation Pipeline

A Nextflow pipeline for detecting and annotating variants from single-cell sequencing data.

## Description

This pipeline processes single-cell sequencing data to identify and annotate genetic variants. It's designed to be reproducible, scalable, and easy to deploy using Docker containers.

## Features

- Quality control of single-cell data
- Variant calling using multiple tools
- Comprehensive variant annotation
- Modular design for easy customization
- Docker support for reproducibility
- Multiple execution profiles (local, docker, cluster, AWS)

## Quick Start

```bash
# Clone the repository
git clone https://github.com/yourusername/sc-variant-annotation-pipeline.git
cd sc-variant-annotation-pipeline

# Run with Docker
nextflow run main.nf -profile docker --input samplesheet.csv
