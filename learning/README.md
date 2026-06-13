# Learning the BioPython library

## Why I started this learning folder

I created this folder as a personal learning space to explore the BioPython library. My motivation comes from a growing interest in bioinformatics and computational biology, and I wanted a hands-on place to experiment with sequence handling, database access, alignment, and analysis using one of the most widely used Python libraries in the field.

The goal is not just to read documentation, but to practice core BioPython workflows and document the kinds of operations that feel most useful for bioinformatics work.

## What I learned here

This learning folder covers the following core BioPython topics:

1. Seq, SeqRecord, and core operations
   - Creating and manipulating `Seq` objects
   - Using `SeqRecord` to store sequence data with annotations and metadata
   - Basic sequence operations like slicing, transcription, translation, and reverse-complement

2. Reading and writing sequence files
   - Loading sequence data from common formats such as FASTA and GenBank
   - Writing sequence objects back to file formats for reuse
   - Understanding `SeqIO` as the main BioPython I/O interface

3. Accessing NCBI databases
   - Querying NCBI nucleotide and protein records
   - Retrieving remote sequence data using Entrez
   - Saving downloaded records for local analysis

4. Running and parsing BLAST
   - Submitting BLAST searches programmatically
   - Parsing BLAST results to extract alignments and scores
   - Integrating BLAST with sequence workflows for variant and similarity discovery

5. Pairwise alignment and variant detection
   - Performing pairwise sequence alignment with `pairwise2`
   - Comparing alignments to find mismatches, insertions, deletions, and variants
   - Using alignments to interpret sequence differences

## Notes about this folder

- This folder is intentionally focused on learning, so examples may be exploratory and incremental.
- Code files in this folder are meant to capture the steps I used while learning each topic.
- The material here is useful for building a practical foundation in BioPython before moving on to more advanced tools like `Bio.Align`, `Bio.Phylo`, or custom analysis pipelines.

## Recommended next steps

- Continue adding examples for `SeqIO` formats such as `GenBank`, `EMBL`, and `Clustal`.
- Add more remote data examples from `Entrez` and `EUtils`.
- Explore `Bio.Align` and multiple sequence alignment.
- Practice building small scripts that combine sequence parsing, BLAST, and alignment into a mini analysis pipeline.

---

This README is a quick reference for why this learning folder exists and what BioPython topics I worked on here.