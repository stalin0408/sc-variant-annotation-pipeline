FROM ubuntu:22.04

LABEL maintainer="stalindanyshaji@gmail.com"
LABEL description="Single-cell variant annotation pipeline"

# Environment
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=UTC
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    wget \
    unzip \
    bzip2 \
    gzip \
    make \
    cmake \
    python3 \
    python3-pip \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libbz2-dev \
    liblzma-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Upgrade pip + install Python dependencies
COPY pyproject.toml /opt/pipeline/
RUN pip3 install --upgrade pip && \
    pip3 install /opt/pipeline/ && \
    pip3 install pandas numpy matplotlib pysam

# Install FastQC
RUN wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip && \
    chmod +x FastQC/fastqc && \
    mv FastQC /opt/ && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm fastqc_v0.12.1.zip

# Install Samtools & BCFtools
RUN wget -q https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar -xjf samtools-1.17.tar.bz2 && \
    cd samtools-1.17 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && make install && \
    cd .. && rm -rf samtools-1.17* && \
    wget -q https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 && \
    tar -xjf bcftools-1.17.tar.bz2 && \
    cd bcftools-1.17 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && make install && \
    cd .. && rm -rf bcftools-1.17*

# Install R packages for QC
RUN R -e "install.packages(c('ggplot2','dplyr','tidyr','plotly','htmlwidgets'), repos='https://cloud.r-project.org/')"

# Copy pipeline files
COPY bin/ /usr/local/bin/
COPY modules/ /opt/pipeline/modules/
COPY subworkflows/ /opt/pipeline/subworkflows/
COPY conf/ /opt/pipeline/conf/
COPY main.nf /opt/pipeline/
COPY nextflow.config /opt/pipeline/

# Create working directories
RUN mkdir -p /data /results /work

WORKDIR /data

# Remove entrypoint (important for Nextflow)
ENTRYPOINT []
CMD ["bash"]