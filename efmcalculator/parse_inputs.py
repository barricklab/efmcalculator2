import pathlib
import logging
import csv
import hashlib
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
from efmcalculator import features
from .constants import FASTA_EXTS, GBK_EXTS
from .EFMSequence import EFMSequence
from .bad_state_mitigation import BadSequenceError, detect_special_cases

logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())


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
    with open(filepath) as f:
        fileinfo = f.read().encode()
    for i, seq in enumerate(sequences):
        try:
            if use_filename:
                filename = Path(filepath).stem
                if not seq.name:
                    seq.name = f"{filename}"
                if not seq.description or seq.description == '':
                    seq.description = f"{filename}"
            originhash = hashlib.md5(fileinfo + bytes(i)).hexdigest()
            efmseq = EFMSequence(seq, originhash)
            test_sequences.append(efmseq)
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
    sequence = str(seq.seq)
    seq_set = set(sequence.upper().replace(" ", ""))
    if not seq_set.issubset(IUPAC_BASES):
        raise BadSequenceError(
            f"Input sequence contains invalid characters. Only IUPAC bases are allowed."
        )
    if "t" in seq_set and "u" in seq_set:
        raise BadSequenceError(
            f"Input sequence has both T and U. This is probably a mistake."
        )
    if seq_set.issubset("RYSWKMBDHVN"):
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
