# DNA to RNA to Protein Translator

This Python program simulates the central dogma of molecular biology: transcription of DNA to RNA, and translation of RNA to proteins. It uses class-based design and file input to process and translate DNA sequences into corresponding protein strings.

## Features

- `DNA` class:
  - Validates DNA sequences
  - Converts DNA to RNA
  - Finds reverse complements

- `RNA` class:
  - Validates RNA sequences
  - Translates RNA into protein sequences using codons

- Supports reading from files where each line contains a frame number and a DNA sequence.

## Example Input File

Each line should include a reading frame and a DNA sequence, like so:

0 ATGCGTACGTTA
1 TGCAGTATGCCG


## How to Run

Run the script with a properly formatted `DNASequences.txt` file:

```bash
python dna_rna_translator.py
