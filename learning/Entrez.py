"""
##Accessing NCBI databases with Biopython's Entrez module##
Entrez gives programmatic access to all NCBI databases
The workflow is always: esearch() to get IDs, then efetch() to download records
Always set your email as a parameter to Entrez.email
"""

from Bio import Entrez, SeqIO
import time

Entrez.email = "adhshreya00@gmail.com"  #required as NCBI blocks without it
#optional- Entrz.api_key = "your_api_key"  :  raises the rate limit from 3 to 10 requests per second

# 1. SEARCH: get database IDs 
handle = Entrez.esearch(
    db="nucleotide",
    term="Mycobacterium tuberculosis rpoB complete cds",
    retmax=20,        #max IDs to return (default 20, hard max 10000)
    sort="relevance"  #or "pub+date", "length"
)
result = Entrez.read(handle)
handle.close()

ids = result["IdList"]
print(f"Total hits: {result['Count']}, fetched: {len(ids)} IDs")

#output obtained when putting in personal email address: 
#   Total hits: 14281, fetched: 20 IDs


# 2. FETCH — download records 
handle = Entrez.efetch(
    db="nucleotide",
    id=",".join(ids[:5]),  #fetch up to 5 at once (comma-separated)
    rettype="fasta",       #fasta | gb | genbank | xml
    retmode="text"
)
records = list(SeqIO.parse(handle, "fasta"))
handle.close()

#Fetch as GenBank for feature annotations
handle = Entrez.efetch(db="nucleotide", id=ids[0], rettype="gb", retmode="text")
gb_record = SeqIO.read(handle, "genbank")
handle.close()


# 3. FETCH PROTEINS
handle = Entrez.esearch(db="protein", term="spike[Gene] SARS-CoV-2", retmax=10)
prot_ids = Entrez.read(handle)["IdList"]
handle.close()

handle = Entrez.efetch(
    db="protein", id=",".join(prot_ids), rettype="fasta", retmode="text"
)
proteins = list(SeqIO.parse(handle, "fasta"))
handle.close()


# 4. DIRECT ACCESSION FETCH 
handle = Entrez.efetch(db="nucleotide", id="NC_045512.2", rettype="gb", retmode="text")
sarscov2 = SeqIO.read(handle, "genbank")
handle.close()
print(f"Fetched: {sarscov2.id}, {len(sarscov2.seq)} bp")


# 5. LOOP SAFELY (respect rate limits) 
for seq_id in ids[:10]:
    handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta", retmode="text")
    rec = SeqIO.read(handle, "fasta")
    handle.close()
    print(f"{rec.id}: {len(rec.seq)} bp")
    time.sleep(0.4)   #stay under 3 req/sec (no API key)

#Databases: nucleotide, protein, pubmed, gene, taxonomy,
#   structure, snp, clinvar, dbvar, sra
