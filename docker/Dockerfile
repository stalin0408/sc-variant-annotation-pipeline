FROM ubuntu:22.04

LABEL maintainer="stalindanyshaji@gmail.com"
LABEL description="Single-cell variant annotation pipeline"

# Set environment variables
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
    && rm -rf /var/lib/apt/lists/*

# Install Python packages
COPY pyproject.toml /opt/pipeline/
RUN pip3 install /opt/pipeline/

# Copy pipeline scripts
COPY bin/ /usr/local/bin/
COPY modules/ /opt/pipeline/modules/
COPY subworkflows/ /opt/pipeline/subworkflows/
COPY main.nf /opt/pipeline/
COPY nextflow.config /opt/pipeline/

WORKDIR /data

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["--help"]
