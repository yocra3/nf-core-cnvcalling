FROM nfcore/base:1.12.1
LABEL authors="Carlos Ruiz-Arenas" \
      description="Docker image containing all software requirements for the nf-core/cnvcalling pipeline"

# Install mamba to manage conda
RUN conda install mamba -c conda-forge
# Install the conda environment
COPY environment.yml /
RUN mamba env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-cnvcalling/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-cnvcalling > nf-core-cnvcalling.yml

## TCAG pipeline
RUN cd ~ && git clone https://github.com/bjtrost/TCAG-WGS-CNV-workflow.git

## ERDS
RUN cd ~ && git clone https://github.com/igm-team/ERDS.git && \
    cd ERDS/erds_tcag/src && make

## MosaicHunter
RUN cd ~ && git clone https://github.com/zzhang526/MosaicHunter.git

# setup aliases
ENV PATH $PATH:/root/ERDS/erds_tcag/src/
