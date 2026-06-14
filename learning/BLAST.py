"""
Running and parsing BLAST
NCBIWWW.qblast() submits a remote job to NCBI (30 sec–5 min)
Always save the XML result to disk (re-running is slow and wastes quota)
Parse with NCBIXML; for production pipelines, install BLAST+ locally instead
"""

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

# REMOTE BLAST (30 sec – 5 min) 
query = SeqIO.read("query.fasta", "fasta")

result_handle = NCBIWWW.qblast(
    program="blastn",        # blastn | blastp | blastx | tblastn | tblastx
    database="nt",           # nt | nr | refseq_rna | swissprot | pdbnt
    sequence=str(query.seq),
    hitlist_size=10,
    expect=1e-5,             # e-value cutoff
    megablast=True           # faster for high-similarity sequences
)

# Save XML immediately — avoids re-running BLAST
with open("blast_result.xml", "w") as f:
    f.write(result_handle.read())


# PARSING RESULTS 
with open("blast_result.xml") as f:
    records = list(NCBIXML.parse(f))

blast_rec = records[0]
print(f"Query: {blast_rec.query[:60]}")
print(f"Database: {blast_rec.database}")
print(f"Hits: {len(blast_rec.alignments)}")

for alignment in blast_rec.alignments:
    for hsp in alignment.hsps:
        if hsp.expect < 1e-5:
            pct_id  = hsp.identities / hsp.align_length * 100
            pct_cov = hsp.align_length / blast_rec.query_length * 100

            print(f"\n  Hit: {alignment.title[:65]}")
            print(f"  Score: {hsp.score:.0f}  |  E-value: {hsp.expect:.1e}")
            print(f"  Identity: {pct_id:.1f}%  |  Coverage: {pct_cov:.1f}%")
            print(f"  Query [{hsp.query_start:>5}]: {hsp.query[:50]}")
            print(f"  Match         : {hsp.match[:50]}")
            print(f"  Sbjct [{hsp.sbjct_start:>5}]: {hsp.sbjct[:50]}")

# LOCAL BLAST (install NCBI BLAST+ for production) 
from Bio.Blast.Applications import NcbiblastnCommandline

blastn_cmd = NcbiblastnCommandline(
    query="query.fasta",
    db="local_db",       # makeblastdb -in seqs.fasta -dbtype nucl -out local_db
    evalue=1e-5,
    outfmt=5,            # XML output
    out="local_result.xml"
)
stdout, stderr = blastn_cmd()   # executes the command

