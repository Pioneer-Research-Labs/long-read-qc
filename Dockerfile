FROM continuumio/miniconda3

RUN conda install -y -c bioconda -c conda-forge \
    seqkit=2.8.2 \
    biopython=1.84 \
    fastp=0.23.4 \
    cutadapt=4.9 \
    levenshtein=0.26.0 \
    regex=2024.9.11 \
    snapgene-reader=0.1.21 \
    scipy=1.14.1 \
    pandas=2.2.2 \
    plotly=5.24.1  \
    click=8.1.7 \
    samtools=1.20 \
    minimap2=2.28 \
    bedtools=2.31.1 \
    seaborn=0.13.2 \
    matplotlib=3.9.2 \
    sourmash=4.8.11 \
    sourmash_plugin_branchwater=0.9.8 \
    fastplong=0.2.2





