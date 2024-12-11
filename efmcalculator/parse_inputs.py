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

    return test_sequences


def validate_sequences(sequences, circular=True, max_len=None):
    cumulative_length = 0
    for seq in sequences:
        cumulative_length += len(seq)
        if max_len and cumulative_length > max_len:
            raise BadSequenceError(
                f"Input sequence(s) is too long. Max length is {max_len} bases."
            )
        validate_sequence(seq, circular=circular)

def validate_sequence(seq, circular=True):
    IUPAC_BASES = set("ACGTURYSWKMBDHVNacgturyswkmbdhvn")
    seq_set = set(seq.upper().replace(" ", ""))
    if not seq_set.issubset(IUPAC_BASES):
        raise BadSequenceError(
            f"Input sequence contains invalid characters. Only IUPAC bases are allowed."
        )
    if "t" in seq_set and "u" in seq_set:
        raise BadSequenceError(
            f"Input sequence has both T and U. This is probably a mistake."
        )
    if seq_set.issubset(RYSWKMBDHVN):
        raise BadSequenceError(
            f"EFM Calculator cannot currently handle IUPAC ambiguity codes."
        )
    detect_special_cases(seq, circular=circular)

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

def detect_special_cases(sequence, circular=True):
    def can_tile(s):
        doubled_s = s + s
        modified_doubled = doubled_s[1:-1]
        return s in modified_doubled
    if can_tile(str(sequence).lower()) and circular:
        raise BadSequenceError("Circular sequence is an infinitely long SSR. Did you mean to use a linear strategy?")
    if len(str(sequence).lower()) < 21 and circular:
        raise BadSequenceError("EFM Calculator is limited to circular sequences of at least 21 bases.")
