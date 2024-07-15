import pathlib
import logging
from Bio import SeqIO

from .constants import FASTA_EXTS, GBK_EXTS

logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())


class BadSequenceError(ValueError):
    pass


def parse_file(filepath: pathlib.Path):
    path_as_string = str(filepath)
    if not filepath.exists():
        raise OSError("File {} does not exist.".format(path_as_string))
    elif filepath.suffix in FASTA_EXTS:
        sequences = SeqIO.parse(path_as_string, "fasta")
    elif filepath.suffix in GBK_EXTS:
        sequences = SeqIO.parse(path_as_string, "genbank")
    else:
        raise ValueError(
            f"File {filepath.suffix} is not a known file format. Must be one of {FASTA_EXTS +GBK_EXTS}."
        )
    print(f"sequences: {sequences}")
    return list(sequences)


def validate_sequences(sequences, max_len=None):
    cumulative_length = 0
    for seq in sequences:
        cumulative_length += len(seq)
        if max_len and cumulative_length > max_len:
            raise BadSequenceError(
                f"Input sequence(s) is too long. Max length is {max_len} bases."
            )
        IUPAC_BASES = set("ACGTURYSWKMBDHVNacgturyswkmbdhvn")
        if not set(seq).issubset(IUPAC_BASES):
            raise BadSequenceError(
                f"Input sequence contains invalid characters. Only IUPAC bases are allowed."
            )
