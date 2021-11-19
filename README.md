# Flu Strain Compare

## Introduction

## Dependencies

Eventually, this will all be containerized, but until then, this requires:

* `python` >= version 3.0
* [`pandas`](https://pypi.org/project/pandas/)
* [`biopython`](https://pypi.org/project/biopython/)
* [`pymol`](https://github.com/schrodinger/pymol-open-source)

## Usage

Once the dependencies are installed, you can modify the `config.json` file depending on what you want to do.

`repo_directory`: Absolute path to the repository root folder.
`figure_directory`: Absolute path to where you want the figures to be saved.
`position_map_infile`: Path to file to convert to standard H3 numbering, relative to the repo root directory.
`seq_file`: Path to fasta-formatted file that contains full-length amino acid HA sequences.
`q1_id`: Sequence ID of the first query strain. The sequence id is the first word in the fasta header of the desired sequence.
`q1_name`: Name of the first query strain.
`q2_id` and `q2_name`: Same as above but for the second query strain.
`seq_lineage`: Just H3N2 for now.

Then, you can run the script by changing into the `src` directory in your terminal and running:

```
pymol make_comparison_figure.py
```

### Strains available
* 2021-2022 southern hemisphere H3N2 cell-based recommendation (name = `A/Darwin/6/2021`, id = `EPI1885402`)
* 2021-2022 northern hemisphere H3N2 Flublok (name = `A/Tasmania/503/2020`, id = `EPI1752480`)
* 2021-2022 northern hemisphere H3N2 cell-based recommendation (name = `A/Cambodia/e0826360/2020`, id =  `EPI1843589`)
* 2020-2021 northern hemisphere H3N2 Flublok (name = `A/Minnesota/41/2019`, id = `EPI1548699`)
* 2020-2021 northern hemisphere H3N2 cell-based recommendation (name = `A/Hong Kong/45/2019`, id = `EPI1409001`) 