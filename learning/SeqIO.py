"""SeqIO handles all file I/O. Use parse() for multi-record files, read() for exactly one
Both return SeqRecord objects
GenBank format is richest as it includes full feature tables, CDS translations, and taxonomy
"""

from Bio import SeqIO

#----------------------READING----------------------#
#Single record (raises ValueError if 0 or 2+ records in file)
record = SeqIO.read("sequence.fasta", "fasta")
print(f"{record.id}: {len(record.seq)} bp")

"""for example, the fasta file could be:
>seq1
ACGTACGTACGT  

the output would be:
seq1: 12 bp
"""

#multiple records: lazy iterator (memory-efficient)
for record in SeqIO.parse("sequences.fasta", "fasta"):
    print(f"{record.id}: {len(record.seq)} bp")

"""for example, the fasta file could be:
>seq1
ACGTACGTACGT
>seq2
GGGCCC

the output would be:
seq1: 12 bp
seq2: 6 bp
"""

#convert to list only when you need random access / small file
records = list(SeqIO.parse("sequences.fasta", "fasta"))

"""for example, the fasta file could be:
>seq1
ACGTACGTACGT
>seq2
GGGCCC

the list would be:
[SeqRecord(seq=Seq('ACGTACGTACGT'), id='seq1', name='seq1', description='seq1', dbxrefs=[]), SeqRecord(seq=Seq('GGGCCC'), id='seq2', name='seq2', description='seq2', dbxrefs=[])]
"""
#GenBank (richest format): features, CDS translations, taxonomy
gb = SeqIO.read("NC_000962.gb", "genbank")
print(gb.annotations["organism"])     #Mycobacterium tuberculosis H37Rv
print(gb.annotations["taxonomy"])     #['Bacteria', 'Actinobacteria', ...]
print(len(gb.features))               #number of annotated features

#extract gene features from GenBank record
for feat in gb.features:
    if feat.type == "CDS":
        gene  = feat.qualifiers.get("gene", ["unknown"])[0]
        prot  = feat.qualifiers.get("translation", [""])[0]
        g_seq = feat.extract(gb.seq)   # handles strand direction automatically
        print(f"{gene}: {len(g_seq)} bp | {prot[:20]}...")

"""for example, the output could be:
rpoB: 3519 bp | MKKLLVLLTALALAGALAAP...
katG: 2223 bp | MKKLLVLLTALALAGALAAP...
"""

#FASTQ: sequences with per-base quality scores
for rec in SeqIO.parse("reads.fastq", "fastq"):
    scores = rec.letter_annotations["phred_quality"]   #list[int]
    avg_q  = sum(scores) / len(scores)
    if avg_q >= 30:
        print(f"{rec.id}: avg Q={avg_q:.1f}")          #keep high-quality reads

"""for example, the output could be:
read1: avg Q=35.2
read2: avg Q=28.7
"""

#----------------------WRITING----------------------#
SeqIO.write(records, "output.fasta", "fasta")

#write a filtered subset
long_seqs = [r for r in records if len(r.seq) > 500]
SeqIO.write(long_seqs, "long.fasta", "fasta")

#format conversion: FASTQ -> FASTA (quality scores are dropped)
SeqIO.write(
    SeqIO.parse("reads.fastq", "fastq"),
    "reads.fasta",
    "fasta"
)

#Supported formats:
# fasta, fastq, genbank (gb), embl, clustal, nexus,
# phylip, stockholm, tab, ig, etc. See Biopython docs for full list