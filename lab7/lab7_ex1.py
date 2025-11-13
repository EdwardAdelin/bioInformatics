#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# lab7_ex1.py
# Detect tandem repeats (unit size 3..6, min repeats 2) in a DNA sequence

import argparse
import random
import sys
import urllib.request
import textwrap
from typing import List, Dict, Tuple


def fetch_fasta_from_ncbi(accession: str) -> str:
    """Fetch a FASTA record from NCBI nuccore using efetch (returns raw fasta text).
    Uses HTTP GET to the eutils efetch endpoint.
    """
    url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}"
        "&rettype=fasta&retmode=text"
    )
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = resp.read().decode("utf-8")
            return data
    except Exception as e:
        raise RuntimeError(f"Failed to fetch accession {accession}: {e}")


def parse_fasta(text: str) -> Tuple[str, str]:
    """Return (header, sequence) from FASTA text. Sequence uppercased and non-ACGTN letters removed."""
    lines = [l.strip() for l in text.splitlines() if l.strip()]
    if not lines:
        return ("", "")
    header = lines[0] if lines[0].startswith(">") else ">unknown"
    seq = "".join(lines[1:] if lines[0].startswith(">") else lines)
    seq = seq.upper()
    # keep only standard letters
    seq = ''.join([c for c in seq if c in "ACGTN"])
    return header, seq


def detect_tandem_repeats(seq: str, min_unit=3, max_unit=6, min_repeats=2) -> List[Dict]:
    """Detect tandem repeats (contiguous repeats) in sequence.
    Returns list of dicts with keys: unit, unit_len, start (0-based), end (exclusive), repeats
    """
    n = len(seq)
    results = []
    for L in range(min_unit, max_unit + 1):
        i = 0
        while i <= n - L * min_repeats:
            unit = seq[i:i+L]
            if len(unit) < L:
                break
            # count how many times unit repeats consecutively
            count = 1
            j = i + L
            while j + L <= n and seq[j:j+L] == unit:
                count += 1
                j += L
            if count >= min_repeats:
                results.append({
                    'unit': unit,
                    'unit_len': L,
                    'start': i,
                    'end': j,
                    'repeats': count,
                })
                # advance to end of this run to avoid reporting sub-runs starting inside
                i = j
            else:
                i += 1
    # sort by start
    results.sort(key=lambda r: (r['start'], -r['unit_len']))
    return results


def format_results(results: List[Dict], seq_name: str, window_info: str = "") -> str:
    if not results:
        return "No tandem repeats (3-6 bp, >=2 repeats) found."
    out_lines = [f"Tandem repeats found in {seq_name} {window_info}:\n"]
    for r in results:
        out_lines.append(
            f"unit={r['unit']} (len={r['unit_len']}), repeats={r['repeats']}, "
            f"positions=[{r['start']}..{r['end']-1}]"
        )
    return "\n".join(out_lines)


def read_fasta_file(path: str) -> Tuple[str, str]:
    with open(path, 'r') as fh:
        text = fh.read()
    return parse_fasta(text)


def pick_random_window(seq: str, min_len=1000, max_len=3000) -> Tuple[int, int, str]:
    n = len(seq)
    if n <= max_len and n >= min_len:
        return 0, n, seq
    if n < min_len:
        raise ValueError(f"Sequence too short ({n} nt) to extract required window of >={min_len} nt")
    L = random.randint(min_len, min(max_len, n))
    start = random.randint(0, n - L)
    return start, start + L, seq[start:start+L]


def main():
    p = argparse.ArgumentParser(description="Detect tandem repeats (3-6 bp) in a DNA sequence window of 1000-3000 nt")
    p.add_argument("--accession", help="NCBI accession to fetch (efetch). If sequence length outside 1000-3000, a random window will be used.")
    p.add_argument("--file", help="Local FASTA file to read")
    p.add_argument("--min-unit", type=int, default=3, help="Minimum repeat unit length (default 3)")
    p.add_argument("--max-unit", type=int, default=6, help="Maximum repeat unit length (default 6)")
    p.add_argument("--min-repeats", type=int, default=2, help="Minimum consecutive repeats (default 2)")
    args = p.parse_args()

    seq_name = "(inline sample)"
    seq = ""
    window_info = ""
    try:
        if args.accession:
            print(f"Fetching accession {args.accession} from NCBI...", file=sys.stderr)
            fasta = fetch_fasta_from_ncbi(args.accession)
            header, full_seq = parse_fasta(fasta)
            seq_name = header
            try:
                s0, s1, seq = pick_random_window(full_seq)
                window_info = f"(window {s0}-{s1}, length={len(seq)})"
                print(f"Using window {s0}-{s1} (len {len(seq)}) from accession {args.accession}", file=sys.stderr)
            except ValueError:
                # full_seq length between 1000-3000
                seq = full_seq
                window_info = f"(full length={len(seq)})"
        elif args.file:
            header, full_seq = read_fasta_file(args.file)
            seq_name = header
            try:
                s0, s1, seq = pick_random_window(full_seq)
                window_info = f"(window {s0}-{s1}, length={len(seq)})"
            except ValueError:
                seq = full_seq
                window_info = f"(full length={len(seq)})"
        else:
            # fallback sample: small synthetic sequence built to be between 1000 and 3000 nt
            # This sequence was generated for example purposes. To use a real NCBI sequence, pass --accession ACCESSION
            seq_name = "sample_synthetic"
            # build a synthetic sequence ~1500 nt with some tandem repeats embedded
            parts = []
            parts.append('A' * 500)
            parts.append(('ATG' * 5))  # a 3-bp unit repeated 5 times
            parts.append('C' * 200)
            parts.append(('GCTA' * 4))  # a 4-bp unit repeated 4 times
            parts.append('T' * 300)
            parts.append(('CAGTAC' * 6))  # a 6-bp unit repeated 6 times
            seq = ''.join(parts)
            window_info = f"(synthetic length={len(seq)})"

        if not seq:
            print("No sequence was loaded.", file=sys.stderr)
            sys.exit(1)

        if len(seq) < 1000 or len(seq) > 3000:
            print(f"Warning: selected sequence length is {len(seq)} (expected 1000-3000).", file=sys.stderr)

        results = detect_tandem_repeats(seq, min_unit=args.min_unit, max_unit=args.max_unit, min_repeats=args.min_repeats)
        out = format_results(results, seq_name, window_info)
        print(out)

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(2)


if __name__ == '__main__':
    main()
