FROM ubuntu:24.04

COPY conda.yml ./

COPY requirements.txt ./

RUN apt-get update && apt-get -y install curl

RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda3-py310_23.10.0-1-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local/ \
    && rm -rf /tmp/miniconda.sh

RUN pip install --no-cache-dir -r requirements.txt

RUN conda-env update -n root -f conda.yml && conda clean -a