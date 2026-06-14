from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

#create seq object 
dna = Seq('ATGCGATCGATCGAATTCG') #dna seq acts as an immutable string

print(dna)  #ATGCGATCGATCGAATTCG
print(len(dna)) #19
print(dna[0:3])  #ATG: start codon

"""Why is start codon inportant? 
The start codon is important because it signals the beginning of translation, which is the process by which the genetic code in mRNA is used to synthesize proteins
The start codon (usually AUG in eukaryotes) indicates where the ribosome should begin translating the mRNA into a protein
It also specifies the amino acid methionine, which is typically the first amino acid incorporated into a newly synthesized protein
Without a start codon, the ribosome would not know where to begin translation, and the resulting protein would not be produced correctly
"""

print(dna.count('AT')) #counts the number of AT in the sequence
print(dna.find('GAATTCG')) #finds the position of the first occurrence of the substring

#biological operations 
print(dna.complement()) #TACGCTAGCTAGCTTAAGC; A->T, T->A, C->G, G->C
print(dna.reverse_complement()) #CGAATTCGATCGATCGCAT; reverse the sequence and then take the complement
print(dna.transcribe()) #AUGCGAUCGAUCGAUUUCG; transcribes DNA to RNA (T->U)
print(dna.translate()) #MRLRSI; translates the DNA sequence to a protein sequence using the standard genetic code

"""Why to translate DNA to protein?
Proteins are the functional molecules in cells that perform a wide variety of tasks, such as catalyzing biochemical reactions, providing structural support, and regulating gene expression
The sequence of amino acids in a protein determines its structure and function, and this sequence is determined by the sequence of nucleotides in the DNA
By translating DNA to protein, we can understand the functional implications of genetic information and how it contributes to the biology of an organism
"""
print(dna.translate(table=11)) #MRLRSI; translates using the bacterial, archaeal and plant plastid code
print(dna.translate(to_stop=True)) #MRLRS; translates until the first stop codon is encountered
#The stop codon is usually represented by an asterisk (*) in the translated protein sequence, and it signals the end of translation. When the to_stop parameter is set to True, the translation process will stop as soon as a stop codon is encountered, and the resulting protein sequence will not include any amino acids that come after the stop codon. This can be useful for analyzing protein sequences that may contain multiple stop codons or for identifying the functional portion of a protein sequence.

#GC content
print(f"{gc_fraction(dna)*100:.1f}%") 

"""why to calculate GC content?
A-T has two bonds, compared to the three bonds between G-C, which is diffciult to break
During certain biological processes, such as DNA replication and transcription, the DNA strands need to be separated
Regions with high GC content are more stable and require more energy to separate, which can affect the efficiency of these processes
Additionally, GC-rich regions can influence the structure of the DNA and its interactions with proteins, impacting gene expression and regulation
Also, GC content can be used to identify and classify organisms, as different species often have characteristic GC content in their genomes
"""

# SeqRecord: sequence + metadata (what SeqIO returns)
record = SeqRecord(
    Seq("ATGCGATCG"),
    id="gene_001",
    name="rpoB",
    description="RNA polymerase beta subunit | M. tuberculosis",
    annotations={"organism": "Mycobacterium tuberculosis",
                 "molecule_type": "DNA"}
)
print(record.id)            # gene_001
print(record.seq)           # ATGCGATCG
print(record.description)   # RNA polymerase beta subunit

# Codon-by-codon iteration
for i in range(0, len(dna) - 2, 3):
    codon = dna[i : i+3]
    if len(codon) == 3:
        print(f"Codon {i//3+1}: {codon} -> {codon.translate()}")
