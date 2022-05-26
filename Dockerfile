FROM ubuntu:20.04
ENV TZ=America/Chicago
RUN apt-get update && \
	ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

RUN apt-get install -y mafft imagemagick libqt5gui5 python3-pip wget
RUN pip install pandas biopython

ENV PYMOL_VERSION 2.5.0

# from https://pymolwiki.org/index.php/Linux_Install#Requirements
RUN apt-get install -y git build-essential python3-dev libglew-dev \
  libpng-dev libfreetype6-dev libxml2-dev \
  libmsgpack-dev python3-pyqt5.qtopengl libglm-dev libnetcdf-dev \
  mafft

RUN wget --no-verbose https://github.com/schrodinger/pymol-open-source/archive/refs/tags/v${PYMOL_VERSION}.tar.gz
RUN tar -xvzf v${PYMOL_VERSION}.tar.gz
WORKDIR pymol-open-source-${PYMOL_VERSION}
RUN python3 setup.py --use-msgpackc=no build install

WORKDIR /app
