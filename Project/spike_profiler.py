"""
SARS-CoV-2 Spike Variant Profiler

A BioPython variant analysis pipeline:
  1. Fetch spike protein sequences from NCBI (Nepal-specific query)
  2. Align each against the Wuhan-Hu-1 reference (NC_045512.2)
  3. Call amino acid mutations via codon-level diffing
  4. Classify sequences against known lineage signatures
  5. Print a formatted report + save to JSON

Requirements:  pip install biopython
Run with:      python spike_profiler.py
"""

from Bio import Entrez, SeqIO, Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from datetime import date
import time, json

Entrez.email = "youremail@example.com"  # required by NCBI

# Variant signature database 
SIGNATURES = {
    "Alpha (B.1.1.7)":   {"N501Y", "P681H", "D614G"},
    "Delta (B.1.617.2)": {"L452R", "T478K", "P681R", "D614G"},
    "Omicron (BA.1)":    {"N501Y", "E484A", "Q493R", "G496S", "D614G", "N679K"},
    "Omicron (BA.2)":    {"N501Y", "T376A", "R408S", "D614G"},
}
REF_ACCESSION = "NC_045512.2"   # Wuhan-Hu-1 reference genome

# 1. Fetch reference spike protein 
def fetch_spike_reference():
    print("[1/4] Fetching reference spike CDS from NCBI...")
    handle = Entrez.efetch(
        db="nucleotide", id=REF_ACCESSION, rettype="gb", retmode="text"
    )
    genome = SeqIO.read(handle, "genbank")
    handle.close()

    for feat in genome.features:
        if feat.type == "CDS":
            product = " ".join(feat.qualifiers.get("product", [])).lower()
            if "spike" in product or "surface glycoprotein" in product:
                spike_nt = feat.extract(genome.seq)
                spike_aa = spike_nt.translate(to_stop=True)
                print(f"    Reference spike: {len(spike_nt)} nt -> {len(spike_aa)} aa")
                return str(spike_aa)
            
    raise ValueError("Spike CDS not found in reference GenBank record")

# 2. Fetch query sequences
def fetch_query_sequences(n=6):
    print(f"[2/4] Fetching spike protein sequences from NCBI...")

    # Nepal-specific query first; fall back to global if no results
    for term in [
        "SARS-CoV-2 spike glycoprotein Nepal",
        "SARS-CoV-2[Organism] spike protein",
    ]:
        handle = Entrez.esearch(db="protein", term=term, retmax=n)
        result = Entrez.read(handle)
        handle.close()
        if result["IdList"]:
            ids = result["IdList"][:n]
            print(f"    Query: '{term}' -> {len(ids)} hits")
            break

    handle = Entrez.efetch(
        db="protein", id=",".join(ids), rettype="fasta", retmode="text"
    )
    seqs = list(SeqIO.parse(handle, "fasta"))
    handle.close()
    print(f"    Downloaded {len(seqs)} sequences")
    return seqs

# 3. Align proteins and call mutations
def find_mutations(ref_aa, query_aa):
    aligner = Align.PairwiseAligner()
    aligner.mode             = "global"
    aligner.open_gap_score   = -10
    aligner.extend_gap_score = -0.5
    try:
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    except Exception:
        aligner.match_score    = 2
        aligner.mismatch_score = -1

    best  = next(iter(aligner.align(ref_aa, query_aa)))
    muts, r_pos = [], 0
    for r, q in zip(str(best[0]), str(best[1])):
        if r != "-": r_pos += 1
        if r != "-" and q != "-" and r != q:
            muts.append(f"{r}{r_pos}{q}")   # e.g. "D614G", "N501Y"
    return muts

# 4. Classify against lineage signatures 
def classify_variant(muts):
    mut_set = set(muts)
    hits = {}
    for lineage, sig in SIGNATURES.items():
        overlap = len(mut_set & sig)
        if overlap:
            hits[lineage] = f"{overlap}/{len(sig)} signature mutations"
    return hits or {"Unknown / Novel lineage": "no signature match"}

# 5. Run full pipeline 
def run_pipeline():
    ref_aa  = fetch_spike_reference()
    queries = fetch_query_sequences(n=6)

    results = []
    print("[3/4] Aligning and calling variants...")

    for seq in queries:
        aa_str = str(seq.seq).replace("-", "").replace("X", "")
        if len(aa_str) < 200:
            continue       # skip fragments too short to be informative

        muts   = find_mutations(ref_aa, aa_str)
        labels = classify_variant(muts)
        results.append({
            "id":          seq.id,
            "description": seq.description[:80],
            "length_aa":   len(aa_str),
            "n_mutations": len(muts),
            "mutations":   muts,
            "lineage":     labels,
        })
        time.sleep(0.3)    # respect NCBI rate limits

    print("[4/4] Generating report...")
    print(f"\n  SARS-CoV-2 SPIKE VARIANT REPORT   {date.today()}\n")
    print(f"  Reference: {REF_ACCESSION}  |  Sequences: {len(results)}")

    for i, r in enumerate(results, 1):
        print(f"  [{i}] {r['id']}")
        print(f"       {r['description']}")
        print(f"       Length: {r['length_aa']} aa  "
              f"|  Mutations vs ref: {r['n_mutations']}")
        if r["mutations"]:
            shown = r["mutations"][:10]
            tail  = f" (+{r['n_mutations']-10} more)" if r["n_mutations"] > 10 else ""
            print(f"       Changes: {', '.join(shown)}{tail}")
        for lineage, detail in r["lineage"].items():
            marker = ">>>" if "Unknown" not in lineage else "   "
            print(f"{marker} Lineage: {lineage} — {detail}")
        print()

    with open("variant_report.json", "w") as f:
        json.dump(results, f, indent=2)
    print("  Full report saved to variant_report.json")


if __name__ == "__main__":
    run_pipeline()