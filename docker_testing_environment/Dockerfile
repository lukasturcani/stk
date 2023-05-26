FROM python:3.11.3-slim

RUN \
  apt-get update \
  && apt-get install -y \
  git \
  make \
  mongo-tools


# requirements.txt is generated with pip-comile. See the justfile.
COPY requirements.txt /deps/requirements.txt

RUN python -m venv /deps/venv && . /deps/venv/bin/activate && pip install -r /deps/requirements.txt
