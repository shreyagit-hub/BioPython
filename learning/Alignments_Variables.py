"""
Bio.Align.PairwiseAligner (BioPython ≥ 1.78) is the modern aligner
After alignment, we iterate over paired characters to find mismatches. These are the SNPs
Codon arithmetic then predicts whether each SNP is synonymous (silent) or nonsynonymous (missense/nonsense)
This is the core of clinical variant analysis
"""

from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq

# PAIRWISE ALIGNMENT (PairwiseAligner — BioPython >= 1.78) 
aligner = Align.PairwiseAligner()
aligner.mode           = "global"  # global: full-length  |  local: best sub-region
aligner.match_score    = 2
aligner.mismatch_score = -1
aligner.open_gap_score   = -5
aligner.extend_gap_score = -0.5

ref = Seq("ATGCGATCGATCGAATCG")
qry = Seq("ATGCGATCGATCGTTCGA")  # 2 SNPs vs reference

best = next(iter(aligner.align(ref, qry)))

print(f"Alignment score: {best.score}")
print(best)                        # pretty-prints the alignment

ref_aligned = str(best[0])         # reference with gap chars ("-")
qry_aligned = str(best[1])         # query with gap chars ("-")

# SNP DETECTION 
def find_snps(ref_seq, query_seq, aligner):
    """Return all SNPs as [{pos (1-indexed), ref, alt}]."""
    best = next(iter(aligner.align(ref_seq, query_seq)))
    snps, ref_pos = [], 0
    for r, q in zip(str(best[0]), str(best[1])):
        if r != "-": ref_pos += 1
        if r != "-" and q != "-" and r != q:
            snps.append({"pos": ref_pos, "ref": r, "alt": q})
    return snps

snps = find_snps(ref, qry, aligner)
for s in snps:
    print(f"SNP @ pos {s['pos']}: {s['ref']} -> {s['alt']}")

# AMINO ACID CHANGE PREDICTION 
def predict_aa_change(snp, ref_nt_seq):
    """Predict synonymous vs missense from SNP position in coding sequence."""
    i            = snp["pos"] - 1         # 0-indexed
    codon_start  = (i // 3) * 3
    codon_num    = (i // 3) + 1           # 1-indexed codon number
    pos_in_codon = i % 3                  # 0, 1, or 2

    ref_codon = str(ref_nt_seq[codon_start : codon_start + 3])
    alt_codon = ref_codon[:pos_in_codon] + snp["alt"] + ref_codon[pos_in_codon+1:]

    ref_aa = str(Seq(ref_codon).translate())
    alt_aa = str(Seq(alt_codon).translate())

    return {
        "codon":      codon_num,
        "notation":   f"{ref_aa}{codon_num}{alt_aa}",  # e.g. "S450L"
        "synonymous": ref_aa == alt_aa,
    }

for snp in snps:
    ch = predict_aa_change(snp, ref)
    tag = "silent" if ch["synonymous"] else "MISSENSE"
    print(f"  {ch['notation']} [{tag}]")

# MULTIPLE SEQUENCE ALIGNMENT 
from Bio import AlignIO

msa = AlignIO.read("alignment.fasta", "fasta")
print(f"MSA: {len(msa)} sequences x {msa.get_alignment_length()} columns")

# Scan columns for polymorphic positions
for col_i in range(msa.get_alignment_length()):
    col      = msa[:, col_i]           # string of one char per sequence
    residues = set(col) - {"-"}
    if len(residues) > 1:
        print(f"Variant at column {col_i+1}: {col}")

# PROTEIN ALIGNMENT — use BLOSUM62 
prot_aligner = Align.PairwiseAligner()
prot_aligner.mode             = "global"
prot_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
prot_aligner.open_gap_score   = -10
prot_aligner.extend_gap_score = -0.5
# then: prot_aligner.align(protein_seq_a, protein_seq_b)
