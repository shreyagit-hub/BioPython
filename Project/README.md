# SARS-CoV-2 Spike Variant Profiler

This project contains a small BioPython-based pipeline for fetching SARS-CoV-2 spike protein sequences from NCBI, comparing them against the Wuhan-Hu-1 reference spike protein, calling amino-acid mutations, and classifying each sequence against known variant signatures.

## Project structure

- `spike_profiler.py` - Main script that runs the pipeline and writes results to `variant_report.json`.
- `variant_report.json` - Generated JSON report containing per-sequence mutation calls and lineage signature matches.

## What it does

The pipeline performs these steps:

1. Fetch the spike coding sequence (CDS) from the Wuhan-Hu-1 reference genome (`NC_045512.2`).
2. Retrieve up to 6 SARS-CoV-2 spike sequence records from NCBI using a sample search term.
3. Align each query protein to the reference spike protein.
4. Call amino-acid substitutions from the alignment.
5. Compare called substitutions against known variant lineage signature mutations.
6. Print a summary report and save the full results to `variant_report.json`.

## Requirements

- Python 3.7+
- BioPython

Install the required package with:

```bash
pip install biopython
```

## Usage

From the `Project` directory, run:

```bash
python spike_profiler.py
```

The script uses `Bio.Entrez` and requires a valid email address in `spike_profiler.py` for NCBI access.

## Generated output

Running the script produces `variant_report.json` in the same folder. The JSON file contains an array of variant records. Each record includes:

- `id`: the sequence accession identifier.
- `description`: the FASTA description text (trimmed to 80 characters).
- `length_aa`: the amino-acid length of the aligned sequence.
- `n_mutations`: the number of amino-acid substitutions versus the Wuhan spike reference.
- `mutations`: a list of mutation strings such as `D614G`, `N501Y`, and `P681H`.
- `lineage`: a mapping of known variant names to signature match counts.

### Example record

```json
{
  "id": "XBA21055.1",
  "description": "XBA21055.1 surface glycoprotein [Severe acute respiratory syndrome coronavirus 2",
  "length_aa": 1176,
  "n_mutations": 22,
  "mutations": [
    "A67V",
    "T95I",
    "G142D",
    "L212I",
    "G339D",
    "S371L",
    "S373P",
    "S375F",
    "K417N",
    "T547K",
    "D614G",
    "F643L",
    "H655Y",
    "N679K",
    "P681H",
    "A701V",
    "N764K",
    "D796Y",
    "N856K",
    "Q954H",
    "N969K",
    "L981F"
  ],
  "lineage": {
    "Alpha (B.1.1.7)": "2/3 signature mutations",
    "Delta (B.1.617.2)": "1/4 signature mutations",
    "Omicron (BA.1)": "2/6 signature mutations",
    "Omicron (BA.2)": "1/4 signature mutations"
  }
}
```

## Notes on the current `variant_report.json`

The current generated report includes multiple spike protein records such as `XBA21055.1`, `XBA21043.1`, `XBA21032.1`, `XBA21020.1`, `XBA21008.1`, and `XBA20996.1`. Each record contains identified mutations and lineage signature overlap counts for the defined variants:

- Alpha (B.1.1.7)
- Delta (B.1.617.2)
- Omicron (BA.1)
- Omicron (BA.2)

When a sequence does not share any configured lineage signatures, the script reports it as `Unknown / Novel lineage`.

## How to extend

- Add more variant signatures to the `SIGNATURES` dictionary in `spike_profiler.py`.
- Adjust the NCBI search query or `retmax` value in `fetch_query_sequences()`.
- Improve mutation calling using full codon-aware alignment or specialized variant calling logic.
