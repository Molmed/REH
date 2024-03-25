# Docker inheritance
FROM bioconductor/bioconductor_docker:devel

# Install required Bioconductor package
ADD ./requirements.txt ./
RUN R -e 'package_list <- read.table("requirements.txt", colClasses="character")[[1]]; \
        BiocManager::install(package_list); \
        devtools::install_github("VanLoo-lab/ascat/ASCAT")'

ENV PATH /miniconda3/bin:$PATH

#    /miniconda3/bin/conda init bash && \
RUN mkdir -p /miniconda3 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda3/miniconda.sh && \
    bash /miniconda3/miniconda.sh -b -u -p /miniconda3 && \
    rm -rf /miniconda3/miniconda.sh && \
    /miniconda3/bin/conda config --add channels defaults && \
    /miniconda3/bin/conda config --add channels bioconda && \
    /miniconda3/bin/conda config --add channels conda-forge && \
    /miniconda3/bin/conda config --set channel_priority strict && \
    /miniconda3/bin/conda install -c bioconda samtools && \
    /miniconda3/bin/conda install -c bioconda tabix && \
    /miniconda3/bin/conda install -c bioconda cancerit-allelecount && \
    echo ". /miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc
