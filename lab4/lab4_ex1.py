"""lab4_ex1.py

Find the coding region of a gene (first ATG .. first in-frame stop codon)
in a DNA sequence, convert the coding region to RNA (T -> U) and translate
it into an amino-acid sequence using the standard genetic code.

Usage examples:
  python lab4_ex1.py ATGAAATTTGGGTAATAG
  python lab4_ex1.py --file example.fasta
  echo ">seq\nATGAAATTTGGGTAATAG" > tmp.fasta; python lab4_ex1.py tmp.fasta

The script accepts either a raw sequence string or a path to a FASTA/text
file. FASTA headers (lines starting with ">") are ignored when reading files.

Behavior:
- Finds the first start codon (ATG). From there, searches for the first
  in-frame stop codon (TAA, TAG, TGA). If a stop codon is found, the coding
  region includes both start and stop codons. If no in-frame stop is found,
  the sequence from the first ATG to the end is used.
- Converts the coding DNA to RNA (T->U) and translates codons to single-letter
  amino-acid codes. Stop codon is represented as '*' at the end of the protein
  sequence if present.

This is a minimal, self-contained implementation intended for lab usage.
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional, Tuple

# Standard genetic code (RNA codons -> single-letter amino acids)
GENETIC_CODE = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}


def clean_sequence(s: str) -> str:
	"""Remove whitespace and FASTA headers and return uppercase DNA sequence."""
	lines = [ln.strip() for ln in s.splitlines()]
	seq_parts = []
	for ln in lines:
		if not ln:
			continue
		if ln.startswith('>'):
			continue
		seq_parts.append(ln)
	seq = ''.join(seq_parts)
	return seq.upper()


def read_input(path_or_seq: Optional[str]) -> str:
	"""If path_or_seq is a path to an existing file, read it. Otherwise
	treat it as a raw sequence. If None, read from stdin."""
	if path_or_seq is None:
		data = sys.stdin.read()
		return clean_sequence(data)

	# try to read file
	try:
		with open(path_or_seq, 'r', encoding='utf8') as fh:
			data = fh.read()
			return clean_sequence(data)
	except FileNotFoundError:
		# treat as raw sequence
		return clean_sequence(path_or_seq)


def find_coding_region(dna: str) -> Optional[Tuple[int, int]]:
	"""Find first ATG start and first in-frame stop (TAA,TAG,TGA).

	Returns (start_idx, end_idx) where indices are 0-based and end_idx is
	inclusive (index of the last base of the stop codon). If no start found,
	returns None. If start found but no in-frame stop found, returns (start, len(dna)-1).
	"""
	dna = dna.upper()
	start = dna.find('ATG')
	if start == -1:
		return None

	stop_codons = {'TAA', 'TAG', 'TGA'}
	# scan in-frame
	for i in range(start + 3, len(dna) - 2, 3):
		codon = dna[i:i+3]
		if codon in stop_codons:
			return start, i+2

	# no in-frame stop found
	return start, len(dna) - 1


def dna_to_rna(dna: str) -> str:
	return dna.replace('T', 'U')


def translate_rna(rna: str) -> str:
	"""Translate an RNA sequence (assumed to start at coding frame).

	Stops are translated to '*' and translation stops at the first stop codon
	if present (including the '*' in the output). If the last codon is
	incomplete it is ignored.
	"""
	prot = []
	for i in range(0, len(rna) - 2, 3):
		codon = rna[i:i+3]
		aa = GENETIC_CODE.get(codon, 'X')
		prot.append(aa)
		if aa == '*':
			break
	return ''.join(prot)


def main(argv: Optional[list] = None) -> int:
	p = argparse.ArgumentParser(description='Translate coding region (ATG..stop) from DNA to protein')
	p.add_argument('input', nargs='?', help='Sequence string or path to a FASTA/text file. If omitted reads from stdin')
	p.add_argument('--no-stop-asterisk', dest='no_stop', action='store_true',
				   help='Do not include "*" for stop codon in protein output')
	args = p.parse_args(argv)

	dna = read_input(args.input)
	if not dna:
		print('No sequence input (empty).', file=sys.stderr)
		return 2

	region = find_coding_region(dna)
	if region is None:
		print('No start codon (ATG) found in input sequence.', file=sys.stderr)
		return 3

	start, end = region
	coding_dna = dna[start:end+1]
	coding_rna = dna_to_rna(coding_dna)
	protein = translate_rna(coding_rna)
	if args.no_stop and protein.endswith('*'):
		protein = protein[:-1]

	# Print results
	print(f'>coding_region {start+1}-{end+1} (1-based indices)')
	print(coding_dna)
	print('\n# RNA (coding)')
	print(coding_rna)
	print('\n# Protein')
	print(protein)
	return 0


if __name__ == '__main__':
	raise SystemExit(main())

