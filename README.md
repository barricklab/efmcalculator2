[![Status](https://github.com/barricklab/efm-calculator2/actions/workflows/package_and_test.yml/badge.svg)](https://github.com/barricklab/efm-calculator2/actions/workflows/package_and_test.yml)

`efmcalculator` is a Python package or web tool for detecting mutational hotspots. It predicts the mutation rates associated with each hotspot and combines them into a relative instability score. These hotspots include simple sequence repeats, repeat mediated deletions, and short repeat sequences. This code updates and improves upon the last version of the [EFM calculator](https://github.com/barricklab/efm-calculator).

`efmcalculator` supports multifasta, genbank, or csv files as input and accepts parameters from the command line. It also supports the scanning of both linear and circular sequences. It defaults to a pairwise comparison strategy (all occurrences of a repeat are compared with all other occurrences), but it also contains an option for a linear comparison strategy (each occurrence of a repeat is only compared with the next occurrence in the sequence) to accelerate the analysis of large sequences.


# Installation
The EFM Calculator can be accessed as a free web tool at efm2-beta.streamlit.app. It is limited to 50000 bases to ensure the app remains performant for other users.
It can be installed and run locally below without such base restriction.

## From pip:
`pip install efmcalculator` or clone this repository and `pip install ./` from the root of the repository.

# Command Line Usage
- -h: help
- -i: inpath
- -o: outpath
- -s: strategy. Either “linear” or “pairwise”
- -c: circular inputs
- -f: output filetype for tables, either csv or parquet
- -j: threads
- -t: tall. Parallelizes across inputs rather than within.
- -v: verbose. 0 (silent), 1 (basic information), 2 (debug)
- --summary: saves only aggrigate results, useful for very tall inputs

Print efmcalculator help:
```
efmcalculator -h
```

Run efmcalculator on all sequences in a FASTA file using the pairwise strategy and print output to csv files within an output folder:
```
efmcalculator -i “input.fasta” -o “output_folder”
```

Run efmcalculator on all sequences in a FASTA file, outputing to the folder output_folder, while treating the input as circular, searching with a linear pattern, and printing debug information:
```
efmcalculator -i “input.fasta” -o “output_folder” -c -s “linear” -v 2
```
