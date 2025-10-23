#!/usr/bin/env python3
"""Compare codon frequencies between SARS-CoV-2 and Influenza (concatenated segments).

Usage: python codon_compare.py

This script expects:
- lab4/NC_045512.2.fasta (SARS-CoV-2) [already present]
- lab4/influenza_raw.fasta (Influenza segments fetched from NCBI)

It will produce PNG charts in the lab4/ directory and print the top-3 amino acids per genome
and an AI prompt that asks which foods contain less of those amino acids.
"""
from collections import Counter, defaultdict
import os
import sys
import matplotlib.pyplot as plt


CODON_TABLE = {
    # DNA codon table (T instead of U)
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}


def read_fasta(path):
    """Return list of (header, sequence) tuples from a FASTA file."""
    records = []
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    with open(path, 'r') as fh:
        header = None
        seq_lines = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line.upper())
        if header is not None:
            records.append((header, ''.join(seq_lines)))
    return records


def concat_sequences(records):
    return ''.join(seq for (_h, seq) in records)


def codon_counts_from_seq(seq, frame=0):
    seq = seq.replace('\n','').replace(' ','').upper()
    # Use T for DNA. Filter only A,C,G,T.
    valid = set('ACGT')
    counts = Counter()
    for i in range(frame, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            continue
        if set(codon) <= valid:
            counts[codon] += 1
    return counts


def aa_counts_from_codon_counts(codon_counts):
    aa_counts = Counter()
    for codon, cnt in codon_counts.items():
        aa = CODON_TABLE.get(codon, 'X')
        aa_counts[aa] += cnt
    return aa_counts


def top_n(counter, n=10):
    return counter.most_common(n)


def plot_top_codons(counter, title, outpath, n=10):
    top = top_n(counter, n)
    codons = [c for c,_ in top]
    vals = [v for _,v in top]
    plt.figure(figsize=(10,6))
    plt.bar(codons, vals, color='C0')
    plt.title(title)
    plt.xlabel('Codon')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def plot_comparison(covid_counts, flu_counts, outpath, top_k=20):
    # union of top codons from both
    covid_top = [c for c,_ in top_n(covid_counts, top_k)]
    flu_top = [c for c,_ in top_n(flu_counts, top_k)]
    combined = sorted(set(covid_top) | set(flu_top))
    covid_vals = [covid_counts.get(c,0) for c in combined]
    flu_vals = [flu_counts.get(c,0) for c in combined]
    x = range(len(combined))
    width = 0.4
    plt.figure(figsize=(max(10, len(combined)*0.4),6))
    plt.bar([i - width/2 for i in x], covid_vals, width=width, label='COVID-19')
    plt.bar([i + width/2 for i in x], flu_vals, width=width, label='Influenza (concatenated)')
    plt.xticks(x, combined, rotation=90)
    plt.ylabel('Count')
    plt.title('Codon counts comparison (top codons union)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


def main():
    here = os.path.dirname(__file__)
    covid_path = os.path.join(here, 'NC_045512.2.fasta')
    influenza_path = os.path.join(here, 'influenza_raw.fasta')

    if not os.path.exists(covid_path):
        print('Missing SARS-CoV-2 fasta at', covid_path, file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(influenza_path):
        print('Missing Influenza fasta at', influenza_path, file=sys.stderr)
        print('Download the influenza segments to', influenza_path, 'and re-run')
        sys.exit(1)

    covid_records = read_fasta(covid_path)
    covid_seq = concat_sequences(covid_records)
    flu_records = read_fasta(influenza_path)
    flu_seq = concat_sequences(flu_records)

    covid_codon_counts = codon_counts_from_seq(covid_seq, frame=0)
    flu_codon_counts = codon_counts_from_seq(flu_seq, frame=0)

    # Top 10 codons plots
    plot_top_codons(covid_codon_counts, 'Top 10 codons - SARS-CoV-2', os.path.join(here, 'covid_top10_codons.png'), n=10)
    plot_top_codons(flu_codon_counts, 'Top 10 codons - Influenza (concatenated)', os.path.join(here, 'flu_top10_codons.png'), n=10)

    # Comparison plot
    plot_comparison(covid_codon_counts, flu_codon_counts, os.path.join(here, 'top_codons_comparison.png'), top_k=10)

    # Amino acid counts
    covid_aa = aa_counts_from_codon_counts(covid_codon_counts)
    flu_aa = aa_counts_from_codon_counts(flu_codon_counts)

    print('\nTop 3 amino acids - SARS-CoV-2:')
    for aa, cnt in covid_aa.most_common(3):
        print(f'  {aa}\t{cnt}')

    print('\nTop 3 amino acids - Influenza (concatenated):')
    for aa, cnt in flu_aa.most_common(3):
        print(f'  {aa}\t{cnt}')

    # Formulate AI prompts using the top-3 amino acids from each genome
    covid_top3 = [aa for aa,_ in covid_aa.most_common(3)]
    flu_top3 = [aa for aa,_ in flu_aa.most_common(3)]

    covid_prompt = (
        "Which common foods are low in the following amino acids? "
        f"{', '.join(covid_top3)}. Provide serving-size examples and approximate gram amounts where possible."
    )

    flu_prompt = (
        "Which common foods are low in the following amino acids? "
        f"{', '.join(flu_top3)}. Provide serving-size examples and approximate gram amounts where possible."
    )

    print('\nAI prompt for SARS-CoV-2 top-3 amino acids:')
    print(covid_prompt)
    print('\nAI prompt for Influenza top-3 amino acids:')
    print(flu_prompt)

    print('\nCharts saved to:')
    print('  -', os.path.join(here, 'covid_top10_codons.png'))
    print('  -', os.path.join(here, 'flu_top10_codons.png'))
    print('  -', os.path.join(here, 'top_codons_comparison.png'))


if __name__ == '__main__':
    main()
