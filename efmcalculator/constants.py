VALID_STRATEGIES = ["pairwise", "linear"]
MIN_SHORT_SEQ_LEN = 0
MAX_SHORT_SEQ_LEN = 16
MIN_SRS_LEN = 6 # Dont search for SRSs below this length
MAX_SSR_SPACER = 5000 # Dont search for SRSs this far appart
UNKNOWN_REC_TYPE = "unknown"
FASTA_EXTS = [".fa", ".fasta"]
GBK_EXTS = [".gb", ".gbk", ".gbff"]
CSV_EXTS = [".csv"]
VALID_EXTS = FASTA_EXTS + GBK_EXTS + CSV_EXTS
SUB_RATE = float(2.2 * 10 ** (-10))
THRESHOLD = 2.2*10e-10 # Minimum mutation rate for a reported hotspot
MAX_SIZE = 50000
MARKER_HEIGHT = 500
VALID_FILETYPES = ["csv", "parquet"]

COLORS = {
    "ori": "#4e7fff",
    "promoter": "#f6a35e",
    "cds": "#479f71",
    "misc": "#808080",
    "primer_bind": "#d0d0d0",
    "terminator": "#C97064",
    "ncrna": "#E8DAB2",
}
