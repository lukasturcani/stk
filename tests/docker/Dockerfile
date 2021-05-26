FROM ubuntu:20.04

ENV PATH="/root/miniconda3/bin:${PATH}"

RUN \
    apt-get update \
    && apt-get install -y \
    curl \
    git \
    make \
    mongo-tools

RUN \
    curl \
    https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh \
    --output miniconda.sh \
    && bash miniconda.sh -b \
    && rm /miniconda.sh

COPY environment.yml /environment.yml

RUN \
    conda env create --name stk --file /environment.yml \
    && rm /environment.yml
