FROM nfcore/base:1.9
LABEL authors="Maxime Borry" \
    description="Docker image containing bio software requirements for the maxibor/metaphlan-nf pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name metaphlan-nf-bio > metaphlan-nf-bio.yml
ENV PATH /opt/conda/envs/metaphlan-nf-bio/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name metaphlan-nf-bio > metaphlan-nf-bio.yml
