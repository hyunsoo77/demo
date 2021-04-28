# Set the base image to CentOS 7.6.1810
FROM centos:7.6.1810

# File Author / Maintainer
MAINTAINER Hyunsoo Kim

# System packages
RUN yum install -y which less sudo vim curl wget parallel bzip2

# Install miniconda3 to /miniconda
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && rm -f Miniconda3-latest-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}

RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda install -c bioconda samtools -y
RUN conda install -c bioconda bcftools -y
RUN conda install -c conda-forge natsort numpy matplotlib psutil seaborn ujson -y

RUN mkdir -p /home/hkim/in /home/hkim/out

WORKDIR /home/hkim
COPY python_plot.py .
RUN ["chmod", "+x", "python_plot.py"]

ENTRYPOINT ["python","/home/hkim/python_plot.py"]



