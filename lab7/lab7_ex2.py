#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lab7_ex2.py
Download 10 influenza genomes from NCBI and, for each genome, detect tandem
repeats (unit lengths 3..6) and plot the most frequent motifs as a bar chart.

Usage examples:
    python3 lab7_ex2.py               # default search term and 10 genomes
    python3 lab7_ex2.py --n 5 --top 8
    python3 lab7_ex2.py --accessions NC_002017,NC_002018  # provide comma-separated accessions

Charts are saved in the lab7/plots directory.
"""

from __future__ import annotations
import argparse
import json
import os
import random
import sys
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

try:
    import matplotlib.pyplot as plt
except Exception as e:
    print("matplotlib is required to run this script. Install it with: pip install matplotlib", file=sys.stderr)
    raise


NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def esearch(term: str, db: str = "nuccore", retmax: int = 10) -> List[str]:
    params = {
        'db': db,
        'term': term,
        'retmax': str(retmax),
        'retmode': 'json'
    }
    url = f"{NCBI_EUTILS_BASE}/esearch.fcgi?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = resp.read().decode('utf-8')
    j = json.loads(data)
    idlist = j.get('esearchresult', {}).get('idlist', [])
    if not idlist:
        # Print raw response to stderr to aid debugging (network, rate-limit, or query syntax issues)
        print(f"esearch returned no IDs for term: {term}", file=sys.stderr)
        print("esearch raw response:\n", data, file=sys.stderr)
    return idlist


def efetch_fasta_by_id(id_: str) -> str:
    url = f"{NCBI_EUTILS_BASE}/efetch.fcgi?db=nuccore&id={urllib.parse.quote(id_)}&rettype=fasta&retmode=text"
    with urllib.request.urlopen(url, timeout=30) as resp:
        return resp.read().decode('utf-8')


def parse_fasta(text: str) -> Tuple[str, str]:
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if not lines:
        return '', ''
    header = lines[0] if lines[0].startswith('>') else '>unknown'
    seq = ''.join(lines[1:] if lines[0].startswith('>') else lines).upper()
    seq = ''.join([c for c in seq if c in 'ACGTN'])
    return header, seq


def detect_tandem_repeats_runs(seq: str, min_unit: int = 3, max_unit: int = 6, min_repeats: int = 2) -> List[Dict]:
    """Detect contiguous tandem repeat runs and return a list of runs.
    Each run is a dict: {'unit': str, 'unit_len': int, 'start': int, 'end': int, 'repeats': int}
    """
    n = len(seq)
    results = []
    for L in range(min_unit, max_unit + 1):
        i = 0
        while i <= n - L * min_repeats:
            unit = seq[i:i+L]
            if len(unit) < L:
                break
            count = 1
            j = i + L
            while j + L <= n and seq[j:j+L] == unit:
                count += 1
                j += L
            if count >= min_repeats:
                results.append({'unit': unit, 'unit_len': L, 'start': i, 'end': j, 'repeats': count})
                i = j
            else:
                i += 1
    results.sort(key=lambda r: (r['start'], -r['unit_len']))
    return results


def aggregate_motif_counts(runs: List[Dict]) -> Counter:
    """Return a Counter mapping motif->score. Score is number of repeated units (repeats).
    e.g. a run of unit 'ATG' repeated 5 times contributes 5 to Counter['ATG'].
    """
    c = Counter()
    for r in runs:
        c[r['unit']] += r.get('repeats', 1)
    return c


def plot_top_motifs(counter: Counter, top_k: int, title: str, outpath: Path) -> None:
    labels = []
    values = []
    for motif, val in counter.most_common(top_k):
        labels.append(motif)
        values.append(val)
    plt.figure(figsize=(max(6, 0.6 * len(labels)), 4))
    plt.bar(range(len(values)), values, color='tab:blue')
    plt.xticks(range(len(values)), labels, rotation=45, ha='right')
    plt.ylabel('Score (sum of repeats)')
    plt.title(title)
    plt.tight_layout()
    plt.savefig(str(outpath), dpi=150)
    plt.close()


def fetch_and_process_accessions(accessions: List[str], outdir: Path, top_k: int, min_unit=3, max_unit=6, min_repeats=2) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    summary = []
    for acc in accessions:
        try:
            print(f"Fetching {acc}...")
            fasta = efetch_fasta_by_id(acc)
            header, seq = parse_fasta(fasta)
            if not seq:
                print(f"Warning: no sequence for {acc}")
                continue
            runs = detect_tandem_repeats_runs(seq, min_unit=min_unit, max_unit=max_unit, min_repeats=min_repeats)
            counts = aggregate_motif_counts(runs)
            if not counts:
                print(f"No repeats found in {acc} ({header})")
            else:
                title = f"Top {top_k} motifs in {acc}"
                outpath = outdir / f"influenza_{acc}.png"
                plot_top_motifs(counts, top_k, title, outpath)
                print(f"Saved plot for {acc} to {outpath}")
            summary.append((acc, header, len(seq), counts.most_common(top_k)))
        except Exception as e:
            print(f"Error processing {acc}: {e}", file=sys.stderr)
    # print summary
    print('\nSummary:')
    for acc, header, length, top in summary:
        print(f"{acc}: {header} len={length} top={top}")


def main():
    p = argparse.ArgumentParser(description='Download influenza genomes and plot most frequent tandem repeats (3-6 bp).')
    p.add_argument('--n', type=int, default=10, help='Number of genomes to download (default 10)')
    p.add_argument('--term', type=str, default='Influenza A virus[Organism] AND "complete genome"[Title]',
                   help='NCBI search term for esearch')
    p.add_argument('--accessions', type=str, help='Comma separated list of accession IDs to use instead of searching')
    p.add_argument('--outdir', type=str, default='plots', help='Directory to save plots (relative to lab7)')
    p.add_argument('--top', type=int, default=8, help='Top K motifs to plot (default 8)')
    args = p.parse_args()

    # Determine lab7 path (script is inside lab7)
    script_dir = Path(__file__).resolve().parent
    outdir = script_dir / args.outdir

    if args.accessions:
        accessions = [a.strip() for a in args.accessions.split(',') if a.strip()]
    else:
        print(f"Searching NCBI for term: {args.term} (retmax={args.n})...", file=sys.stderr)
        ids = esearch(args.term, retmax=args.n)
        if not ids:
            print('No IDs found from esearch.', file=sys.stderr)
            sys.exit(1)
        accessions = ids[:args.n]

    # fetch and process
    fetch_and_process_accessions(accessions, outdir, args.top)


if __name__ == '__main__':
    main()
