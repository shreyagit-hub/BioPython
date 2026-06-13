"""
Running and parsing BLAST
NCBIWWW.qblast() submits a remote job to NCBI (30 sec–5 min)
Always save the XML result to disk (re-running is slow and wastes quota)
Parse with NCBIXML; for production pipelines, install BLAST+ locally instead
"""

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

