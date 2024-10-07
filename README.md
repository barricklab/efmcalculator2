`efmcalculator` is a Python package or web tool for detecting mutational hotspots. It predicts the mutation rates associated with each hotspot and combines them into a relative instability score. These hotspots include simple sequence repeats, repeat mediated deletions, and short repeat sequences. This code updates and improves upon the last version of the [EFM calculator](https://github.com/barricklab/efm-calculator).

`efmcalculator` supports multifasta, genbank, or csv files as input and accepts parameters from the command line. It also supports the scanning of both linear and circular sequences. It defaults to a pairwise comparison strategy (all occurrences of a repeat are compared with all other occurrences), but it also contains an option for a linear comparison strategy (each occurrence of a repeat is only compared with the next occurrence in the sequence) to accelerate the analysis of large sequences. 


# Installation 
Efmcalculator can be accessed as a free web tool at ____.

## From pip:
Run pip install efmcalculator 

# Command Line Usage
- -h: help
- -i: inpath
- -o: outpath
- -s: strategy. Either “linear” or “pairwise”
- -c: circular
- -v: verbose. 0 (silent), 1 (basic information), 2 (debug)
- --no-vis: skip visualization

Print efmcalculator help:
```
efmcalculator -h
```

Run efmcalculator on all sequences in a FASTA file using the pairwise strategy and print output to csv files within an output folder:
```
efmcalculator -i “input.fasta” -o “output_folder”
```

Run efmcalculator on all circular sequences in a FASTA file using the pairwise strategy and print output to csv files within an output folder:
```
efmcalculator -i “input.fasta” -o “output_folder” -c
```

Run efmcalculator on all sequences in a FASTA file using the linear strategy and print output to csv files within an output folder:
```
efmcalculator -i “input.fasta” -o “output_folder” -s “linear”
```

Run efmcalculator on all sequences in a FASTA file using the pairwise strategy and print output to csv files within an output folder without any visualization:
```
efmcalculator -i “input.fasta” -o “output_folder” --no-vis
```

Run efmcalculator on all sequences in a FASTA file using the pairwise strategy and print output to csv files within an output folder in debugging mode:
```
efmcalculator -i “input.fasta” -o “output_folder” -v 2
```
