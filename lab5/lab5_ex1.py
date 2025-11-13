"""
Enounce:
1. Take an arbitrary DNA sequence from NCBI, between 1000 and 3000 nucleotides.
a) Take 2000 random samples from this sequence, of about 100-150  bases
b) Store these samples in an array variable/list
c) Rebuild the original DNA sequence by using these random samples.
d) What will be your main problem with you algorithm? Describe ; namely what kind of strcuture
inside the original sequence may create different computation issues.
Note: the samples must be aligned starting with the minimum of 10 positions in order to avoid random matching.
"""

"""Implementation notes:
 - This script reads `synthetic_sequence_2000nt.fasta` (in the same folder) if present.
 - It samples 2000 reads of length 100-150 whose start positions are multiples of 10 (to respect the "aligned starting with the minimum of 10 positions" requirement).
 - It stores reads in a list and attempts to reconstruct the original sequence using a greedy overlap assembler that requires a minimum overlap of 10.
 - At the end it prints assembly statistics and a short discussion of problems/edge-cases.
"""

import os
import random
import sys
from collections import defaultdict


def read_fasta(path):
	if not os.path.exists(path):
		return None
	seq = []
	with open(path, "r") as f:
		for line in f:
			if line.startswith(">"):
				continue
			seq.append(line.strip())
	return "".join(seq)


def sample_reads(sequence, n_reads=2000, min_len=100, max_len=150, align=10):
	reads = []
	L = len(sequence)
	starts = [i for i in range(0, L - min_len + 1) if i % align == 0]
	if not starts:
		raise ValueError("No valid start positions with the given alignment and sequence length")
	for _ in range(n_reads):
		start = random.choice(starts)
		rlen = random.randint(min_len, max_len)
		if start + rlen > L:
			# clamp
			start = max(0, L - rlen)
		reads.append(sequence[start:start + rlen])
	return reads


def overlap(a, b, min_length=10):
	"""Return length of longest suffix of 'a' matching prefix of 'b' with at least min_length.
	Checks from longest possible down to min_length.
	"""
	start = min_length
	max_ol = 0
	max_possible = min(len(a), len(b))
	# check decreasing lengths for speed (we want the longest)
	for l in range(max_possible, start - 1, -1):
		if a[-l:] == b[:l]:
			return l
	return 0


def build_prefix_dict(reads, k=10):
	d = defaultdict(set)
	for i, r in enumerate(reads):
		if len(r) >= k:
			d[r[:k]].add(i)
	return d


def greedy_assemble(reads, min_overlap=10):
	"""Greedy assembly using min_overlap. Returns list of contigs after no more merges possible.
	Uses k-mer prefix index of length min_overlap to reduce candidate checks.
	"""
	reads = list(reads)  # copy
	k = min_overlap
	prefix_dict = build_prefix_dict(reads, k)

	while True:
		best_i = None
		best_j = None
		best_ol = 0
		# For each read, only consider candidates whose prefix (k) matches the suffix k of the read
		for i, a in enumerate(reads):
			if a is None:
				continue
			if len(a) < k:
				continue
			suf = a[-k:]
			candidates = prefix_dict.get(suf, set())
			for j in candidates:
				if i == j:
					continue
				b = reads[j]
				if b is None:
					continue
				ol = overlap(a, b, min_length=k)
				if ol > best_ol:
					best_ol = ol
					best_i = i
					best_j = j

		if best_ol >= min_overlap and best_i is not None:
			# merge best_i and best_j
			a = reads[best_i]
			b = reads[best_j]
			merged = a + b[best_ol:]
			# replace a with merged, mark j as removed
			reads[best_i] = merged
			reads[best_j] = None
			# update prefix_dict for merged read
			if len(merged) >= k:
				prefix_dict[merged[:k]].add(best_i)
			# remove any references to best_j from prefix_dict
			for key in list(prefix_dict.keys()):
				if best_j in prefix_dict[key]:
					prefix_dict[key].discard(best_j)
				if not prefix_dict[key]:
					del prefix_dict[key]
			continue
		else:
			break

	# filter out removed reads
	contigs = [r for r in reads if r is not None]
	return contigs


def best_alignment_score(short, long):
    """Slide short along long and return best match count and percent identity.
    Returns (best_matches, percent_identity, best_offset)
    """
    best = 0
    best_off = 0
    for off in range(-len(short) + 1, len(long)):
        matches = 0
        overlap_len = 0
        for i in range(len(short)):
            j = off + i
            if 0 <= j < len(long):
                overlap_len += 1
                if short[i] == long[j]:
                    matches += 1
        if overlap_len > 0 and matches > best:
            best = matches
            best_off = off
    pct = 100.0 * best / max(1, min(len(short), len(long)))
    return best, pct, best_off


# add helper to compute reverse complement
def revcomp(s: str) -> str:
    TRANS = str.maketrans("ACGTacgt", "TGCAtgca")
    return s.translate(TRANS)[::-1]


def main():
    # Try to read the synthetic sequence provided in the lab folder
    here = os.path.dirname(os.path.abspath(__file__))
    fasta_path = os.path.join(here, "synthetic_sequence_2000nt.fasta")
    seq = read_fasta(fasta_path)
    if seq is None:
        print("Warning: no fasta found at", fasta_path)
        print("Please provide a sequence file or update the path. Exiting.")
        return

    # enforce enounce length requirement
    if not (1000 <= len(seq) <= 3000):
        print(f"Error: sequence length {len(seq)} not in required range 1000-3000 nt. Exiting.")
        return

    print("Original sequence length:", len(seq))

    # 1.a/b: sample reads
    random.seed(42)
    reads = sample_reads(seq, n_reads=2000, min_len=100, max_len=150, align=10)
    print("Sampled reads:", len(reads))

    # Optionally save samples
    samples_out = os.path.join(here, "samples_2000.fasta")
    with open(samples_out, "w") as f:
        for i, r in enumerate(reads):
            f.write(f">read_{i}\n")
            f.write(r + "\n")
    print("Wrote samples to:", samples_out)

    # 1.c: Reconstruct
    contigs = greedy_assemble(reads, min_overlap=10)
    contigs_sorted = sorted(contigs, key=len, reverse=True)
    print("Number of contigs after greedy assembly:", len(contigs_sorted))
    print("Longest contig length:", len(contigs_sorted[0]) if contigs_sorted else 0)

    # write contigs to a FASTA so you can inspect the reconstructed sequences
    contigs_out = os.path.join(here, "contigs_greedy_minov10.fasta")
    with open(contigs_out, "w") as cf:
        for i, c in enumerate(contigs_sorted):
            cf.write(f">contig_{i}_len{len(c)}\n")
            cf.write(c + "\n")
    print("Wrote contigs to:", contigs_out)

    # Compare longest contig to original sequence using best alignment (also check reverse complement)
    if contigs_sorted:
        longest = contigs_sorted[0]
        # check both orientations of reference
        bm_fwd = best_alignment_score(longest, seq)
        rc_seq = revcomp(seq)
        bm_rev = best_alignment_score(longest, rc_seq)
        best_matches, pct, off = bm_fwd if bm_fwd[0] >= bm_rev[0] else bm_rev
        print(f"Best alignment matches: {best_matches}; percent identity (relative to shorter) ~ {pct:.2f}%; offset {off}")

    # 1.d: Discuss main problems
    print("\nMain algorithmic problems and sequence structures that cause issues:")
    print("- Repeats: long exact repeats (longer than read length or longer than min overlap) make it ambiguous how to order reads; greedy merges may place repeats wrongly.")
    print("- Low-complexity regions (e.g., homopolymers) create many spurious overlaps and false merges.")
    print("- Insufficient coverage in some regions: if some parts of the original have few or no reads covering them (due to sampling), reconstruction will break into multiple contigs.")
    print("- The greedy algorithm is not optimal: it can make local choices that prevent globally optimal assembly (scaffolding issue).")
    print("- Errors/noise in reads (not present here) would further confuse overlaps and require error-tolerant matching.")


if __name__ == "__main__":
	main()
