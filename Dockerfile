FROM ubuntu:20.04
ENV TZ=America/Chicago
RUN apt-get update && \
	ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone && \
	apt-get install -y libqt5gui5 pymol && \
    apt-get install -y python3-pip &&\
    pip install pandas && \
    pip install biopython
COPY ./src /usr/src
WORKDIR /usr/src
CMD ["pymol", "-c", "./make_comparison_figure.py"]