#FROM continuumio/anaconda3:latest
FROM ensemblorg/ensembl-vep:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

USER root
WORKDIR /root

# Install Dependencies of Miniconda
RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 curl git build-essential libdbi-perl nano && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install miniconda3
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# Set Up Channels
RUN conda update -y conda && \
    conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge

# ENV install
COPY env_install.sh /root/env_install.sh
RUN bash env_install.sh

# vcf2maf install
RUN wget https://github.com/mskcc/vcf2maf/archive/refs/tags/v1.6.21.zip && \
    unzip v1.6.21.zip && \
    rm -rf v1.6.21.zip 

# ANNOVAR install
COPY annovar.latest.tar.gz /root/annovar.latest.tar.gz
RUN tar -xzvf annovar.latest.tar.gz && \
    perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene annovar/humandb/ && \
    perl annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene annovar/humandb/ && \
    perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene annovar/humandb/ && \
    perl annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar ensGene annovar/humandb/ && \
    perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar icgc28 annovar/humandb/ && \
    perl annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar icgc28 annovar/humandb/ && \
    perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150 annovar/humandb/ && \
    perl annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar avsnp150 annovar/humandb/ && \
    perl annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20210501 annovar/humandb/ && \
    perl annovar/annotate_variation.pl -buildver hg38 -downdb -webfrom annovar clinvar_20210501 annovar/humandb/ && \
    rm -rf annovar.latest.tar.gz

CMD [ "/bin/bash" ]