import pathlib
import logging
import csv
from Bio import SeqIO, SeqRecord, Seq
from pathlib import Path

from .constants import FASTA_EXTS, GBK_EXTS

logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())


class BadSequenceError(ValueError):
    pass


def parse_file(filepath: pathlib.Path, use_filename: bool = True) -> list:
    """
    parses the inputted files and returns a list of sequences found in each file

    input:
        filepath: path to the multifasta, genbank, or csv file containing the sequences to be scanned

    returns:
        list containing all the sequences to be scanned
    """

    path_as_string = str(filepath)
    if not filepath.exists():
        raise OSError("File {} does not exist.".format(path_as_string))
    elif filepath.suffix in FASTA_EXTS:
        sequences = SeqIO.parse(path_as_string, "fasta")
    elif filepath.suffix in GBK_EXTS:
        sequences = SeqIO.parse(path_as_string, "genbank")
    elif filepath.suffix in GBK_EXTS:
        sequences = parse_csv(path_as_string)
    else:
        raise ValueError(
            f"File {filepath} is not a known file format. Must be one of {FASTA_EXTS + GBK_EXTS + [".csv"]}."
        )

    # If the genbank file doesnt have a name, add it
    test_sequences = []
    for i, seq in enumerate(sequences):
        try:
            if use_filename:
                filename = Path(filepath).stem
                if not seq.name:
                    seq.name = f"{filename}"
                if not seq.description or seq.description == '':
                    seq.description = f"{filename}"
            test_sequences.append(seq)
        except:
            pass

    for seq in test_sequences:
        print(f"{seq.name=}, {seq.description=}")

    return test_sequences


def validate_sequences(sequences, max_len=None):
    cumulative_length = 0
    for seq in sequences:
        cumulative_length += len(seq)
        if max_len and cumulative_length > max_len:
            raise BadSequenceError(
                f"Input sequence(s) is too long. Max length is {max_len} bases."
            )
        validate_sequence(seq)

def validate_sequence(seq):
    IUPAC_BASES = set("ACGTURYSWKMBDHVNacgturyswkmbdhvn")
    if not set(seq).issubset(IUPAC_BASES):
        raise BadSequenceError(
            f"Input sequence contains invalid characters. Only IUPAC bases are allowed."
        )

def parse_csv(path_as_string):
    with open(path_as_string, "r") as csvfile:
        csvreader = csv.reader(csvfile)
        headers = csvreader.__next__()
        if not headers:
            raise BadSequenceError("CSV file is empty or malformed")
        names = headers.index("name")
        seqs = headers.index("seq")
        if not seqs:
            raise BadSequenceError(
                f"CSV file has no 'csv' column"
            )
        sequences = []
        for i, entry in enumerate(sequences):
            if names:
                name = f"{i+1}_" + entry[names]
            else:
                name = f"{i+1}_sequence"
            sequence = entry[seqs]
            formatted_entry = SeqRecord(Seq(sequence),name=name)
            sequences.append(formatted_entry)
    return sequences
