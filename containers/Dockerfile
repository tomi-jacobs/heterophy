FROM continuumio/miniconda3:23.10.0-1

LABEL maintainer="Your Lab <your@email.edu>"
LABEL description="HeteroPhy: SNP & Heterozygosity-Aware Phylogenetic Inference Pipeline"
LABEL version="1.0.0"

# ── System dependencies ────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
        procps \
        wget \
        curl \
        build-essential \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# ── Conda channels ─────────────────────────────────────────────────────────
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# ── Install all tools in one environment ──────────────────────────────────
COPY envs/ /tmp/envs/

RUN conda create -n heterophy --file /tmp/envs/alignment.yml && \
    conda env update -n heterophy --file /tmp/envs/variants.yml && \
    conda env update -n heterophy --file /tmp/envs/phylo.yml && \
    conda env update -n heterophy --file /tmp/envs/python.yml && \
    conda clean -afy

# ── Add pipeline scripts ───────────────────────────────────────────────────
COPY bin/ /opt/heterophy/bin/
RUN chmod +x /opt/heterophy/bin/*.py && \
    ln -s /opt/heterophy/bin/*.py /usr/local/bin/

ENV PATH="/opt/conda/envs/heterophy/bin:/opt/heterophy/bin:$PATH"

# ── Activate environment by default ───────────────────────────────────────
RUN echo "source activate heterophy" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

WORKDIR /data

CMD ["bash"]
