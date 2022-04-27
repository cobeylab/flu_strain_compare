# Flu Strain Compare

## Purpose

This repo contains utilities to visualize mutations between pairs of HA sequences. `make_comparison_figure.py` takes two HA sequences as input and outputs a figure highlighting amino acid and PNGS changes on a representative HA crystal structure. `make_movie.py` can take an ordered list of HA sequences as input and outputs a movie that depicts the amino acid and PNGS changes that occurred between consecutive pairs of strains in the list. At present, H1pdm and H3 strains are supported.

## Dependencies

This utility runs in a container defined by the `Dockerfile`, so you should only need Docker installed on your system. You can use the following command to build the container, but note that Docker does not run on Midway. I've tested this on Ubuntu, Windows 10, and an M1 Mac.

```
docker build ./ -t flu_strain_compare
```

You can, however, use Singularity on Midway. You have to pull the image from my (Phil's) account on the Singularity cloud since we cannot build de novo containers on the cluster.

```
module load singularity
singularity pull --arch amd64 library://philarevalo/dev/ubuntu-pymol-biopython:latest
```

## Configuration

The `configuration/config.json` file serves as input for the `make_comparison_figure.py` script.

* `seq_file`: Name of fasta-formatted file that contains full-length amino acid HA sequences. File must be in `data` directory.
* `q1_id`: Sequence ID of the first query strain. The sequence id is the first word in the fasta header of the desired sequence.
* `q2_id`: Same as above but for the second query strain.
* `seq_lineage`: Specify the lineage of your query strains. Either H1 or H3 for now.
* `numbering_scheme`: What numbering scheme do you want to use for mutation identification? For H1 sequences, you can choose `H1pdm`, `H3`, or `H1_1933`. For H3 sequences, only `H3` numbering is available.

Similarly, the `configuration/movie_config.json` file serves as input for the `make_movie.py` script. 
* `output_handle`: The output base filename for the final movie.
* `seq_file`: Name fasta-formatted file that contains full-length amino acid HA sequences. File must be in `data` directory.
* `seq_lineage`: Specify the lineage of your query strains. Either H1 or H3 for now.
* `numbering_scheme`: What numbering scheme do you want to use for mutation identification? For H1 sequences, you can choose `H1pdm`, `H3`, or `H1_1933`. For H3 sequences, only `H3` numbering is available.
* `frame_order`: The ordered list of sequence IDs to make the movie. Each frame of the movie consists of a comparison figure between each consecutive pair of strains in the list.

## Running 

With your configuration file set up to your liking, you can just run the container with the following commands:

```
export flu_strain_compare_path=<ABSOLUTE PATH TO REPO DIRECTORY>
docker run -v ${flu_strain_compare_path}:/app flu_strain_compare python3 <SCRIPT NAME>.py
```

For Singularity:

```
export flu_strain_compare_path=<ABSOLUTE PATH TO REPO DIRECTORY>
singularity exec --bind ${flu_strain_compare_path}:/app ubuntu-pymol-biopython_latest.sif pymol -c /app/src/<SCRIPT NAME>.py
```

In either case, `<SCRIPT NAME>` should either be `make_comparison_figure` or `make_movie`/


## Outputs

Outputs from both scripts will be written to the `figures` directory of the repository.

## Strains available
### H3
* 2021-2022 southern hemisphere cell-based recommendation (name = `A/Darwin/6/2021`, id = `EPI1885402`)
* 2021-2022 northern hemisphere Flublok (name = `A/Tasmania/503/2020`, id = `EPI1752480`)
* 2021-2022 northern hemisphere cell-based recommendation (name = `A/Cambodia/e0826360/2020`, id =  `EPI1843589`)
* 2020-2021 northern hemisphere Flublok (name = `A/Minnesota/41/2019`, id = `EPI1548699`)
* 2020-2021 northern hemisphere cell-based recommendation (name = `A/Hong Kong/45/2019`, id = `EPI1409001`) 
* 2019-2020 cell-based recommendation (name = `A/Kansas/14/2017`, id = `EPI1174043`)
* 2018-2019 cell-based recommendation (name = `A/Singapore/INFIMH-16-0019/2016`, id = `EPI1106235`)
* 2016-2018 cell-based recommendation (name = `A/Hong Kong/4801/2014`, id = `EPI539576`)
* 2015-2016 cell-based recommendation (name = `A/Switzerland/9715293/2013`, id = `EPI530687`)
* 2013-2015 cell-based recommendation (name = `A/Victoria/361/2011 `, id = `EPI349103`)
### H1
* 2021-2022 northern hemisphere cell-based recommendation (name = `A/Wisconsin/588/2019`, id = `EPI1715168`)
* 2020-2021 northern hemisphere cell-based recommendation (name = `A/Hawaii/70/2019`, id = `EPI1669665`) 
* 2019-2020 cell-based recommendation (name = `A/Brisbane/02/2018`, id = `EPI1212884`)
* 2017-2019 cell-based recommendation (name = `A/Michigan/45/2015`, id = `EPI699812`
* 2009-2017 cell-based recommendation (name = `A/California/04/2009`, id = `EPI178457`)
* 2013-2014 example from Linderman et al. (name = `A/Colorado/3514/2013`, id = `EPI501723`)
