# Flu Strain Compare

## Introduction

## Dependencies

This utility runs in a container defined by the `Dockerfile`, so you should only need Docker installed on your system to build the container with the following command:

```
docker build ./ -t flu_strain_compare
```

If you're using Singularity, then it's a little different:

```
singularity pull --arch amd64 library://philarevalo/dev/ubuntu-pymol-biopython:latest
```

## Configuration

Once the dependencies are installed, you can modify the `configuration/config.json` file depending on what you want to do.

* `position_map_infile`: Path to file to convert to standard H3 numbering. File must be in `data` directory.
* `seq_file`: Path to fasta-formatted file that contains full-length amino acid HA sequences. File must be in `data` directory.
* `q1_id`: Sequence ID of the first query strain. The sequence id is the first word in the fasta header of the desired sequence.
* `q1_name`: Name of the first query strain.
* `q2_id` and `q2_name`: Same as above but for the second query strain.
* `seq_lineage`: Just H3N2 for now.

## Running 

With your configuration file set up to your liking, running is as simple as

```
export flu_strain_compare_path=<ABSOLUTE PATH TO REPO DIRECTORY>
docker run -v ${flu_strain_compare_path}/figures:/usr/figures -v ${flu_strain_compare_path}/configuration:/usr/configuration -v ${flu_strain_compare_path}/data:/usr/data flu_strain_compare
```

For Singularity:

```
export flu_strain_compare_path=<ABSOLUTE PATH TO REPO DIRECTORY>
singularity exec --bind ${flu_strain_compare_path}/figures:/usr/figures,${flu_strain_compare_path}/configuration:/usr/configuration,${flu_strain_compare_path}/data:/usr/data,${flu_strain_compare_path}/src:/usr/src ubuntu-pymol-biopython_latest.sif pymol -c /usr/src/make_comparison_figure.py
```

### Strains available
* 2021-2022 southern hemisphere H3N2 cell-based recommendation (name = `A/Darwin/6/2021`, id = `EPI1885402`)
* 2021-2022 northern hemisphere H3N2 Flublok (name = `A/Tasmania/503/2020`, id = `EPI1752480`)
* 2021-2022 northern hemisphere H3N2 cell-based recommendation (name = `A/Cambodia/e0826360/2020`, id =  `EPI1843589`)
* 2020-2021 northern hemisphere H3N2 Flublok (name = `A/Minnesota/41/2019`, id = `EPI1548699`)
* 2020-2021 northern hemisphere H3N2 cell-based recommendation (name = `A/Hong Kong/45/2019`, id = `EPI1409001`) 