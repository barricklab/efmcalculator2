VALID_STRATEGIES = ["pairwise", "linear"]
MIN_SHORT_SEQ_LEN = 0
MAX_SHORT_SEQ_LEN = 16
MIN_SSR_LEN = 6 # Dont search for SRSs below this length
MAX_SSR_SPACER = 5000 # Dont search for SRSs this far appart
UNKNOWN_REC_TYPE = "unknown"
FASTA_EXTS = [".fa", ".fasta"]
GBK_EXTS = [".gb", ".gbk", ".gbff"]
VALID_EXTS = FASTA_EXTS + GBK_EXTS
SUB_RATE = float(2.2 * 10 ** (-10))

COLORS = {
    "ori": "#4e7fff",
    "promoter": "#f6a35e",
    "cds": "#479f71",
    "misc": "#808080",
    "primer_bind": "#d0d0d0",
    "terminator": "#C97064",
    "ncrna": "#E8DAB2",
}
