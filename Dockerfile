FROM r-base:4.1.3
ARG DEBIAN_FRONTEND=noninteractive


RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('ggplot2')"
COPY requirements.R requirements.R


RUN apt-get update && \
    apt-get install -y sudo && \
    apt-get install -y jq && \
    apt-get install expect  -y  --no-install-recommends \
        wget \
        zip unzip \
        parallel \
        curl \
        git \
        seqtk \
        && \
    apt-get install -y libcurl4-openssl-dev  \
        default-jre \
        fastqc \
        bedtools \
        samtools \
        procps 


 
RUN Rscript requirements.R

# Install miniconda
ENV CONDA_DIR /opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda 
ENV PATH=$CONDA_DIR/bin:$PATH
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -c bioconda pbtk  

     
RUN apt-get install -y python3.10
RUN apt-get install -y python3-pip && \
    pip install pandas && \
    pip install numpy  && \
    conda install matplotlib && \
    apt-get install -y  procps && \
    pip install pysam && \
    pip install sniffles
RUN pip install configparser && \
    pip install biopython && \
    pip3 install bio --user && \
    apt-get install -y ncbi-blast+ -y

# Install aws-cli
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install
RUN mkdir tools && \
    cd tools && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
    chmod +x blat  && \
    wget -O BBMap_39.01.tar.gz  https://sourceforge.net/projects/bbmap/files/latest/download && \
    tar -xvzf BBMap_39.01.tar.gz && \ 
    wget -O bwakit-0.7.15_x64-linux.tar.bz2  https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2/download && \
    tar -xvjf bwakit-0.7.15_x64-linux.tar.bz2 && \
    wget https://github.com/biod/sambamba/releases/download/v0.8.2/sambamba-0.8.2-linux-amd64-static.gz && \
    gunzip sambamba-0.8.2-linux-amd64-static.gz && \
    mv sambamba-0.8.2-linux-amd64-static sambamba-0.8.2 && \
    chmod +x sambamba-0.8.2 && \
    ln -s sambamba-0.8.2 /usr/local/bin
COPY import import
COPY usearch11.0.667_i86linux64 ./tools/usearch


 
ENTRYPOINT ["bash", "/import/entrypoint.sh"]

