"""
Find the percentage for all the dinucleotide 
and trinucleotide combinations 
for the sequence: 
S="TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA".

1. Build a brute force engine to generate all 
dinucleotide and trinucleotide combinations.

2. For each combination, find out the percentage 
inside the S sequence.

3. Show the percentage for each combination in the 
output of your implementation.
"""

# Exercise nr. 1 and nr. 2
S = "TACGTGCGCGCGAGCTATCTACTGACTTACGACTAGTGTAGCTGCATCATCGATCGA"
from itertools import product # this is used to generate combinations
# the next function generates all possible combinations of kmers
def generate_kmers(k, alphabet=('A', 'C', 'G', 'T')): # k = length, mers = from greek
	return [''.join(p) for p in product(alphabet, repeat=k)] # all possible combinations (length k)

def count_kmers(sequence, k): # counts overlapping kmers and returns a dictionary
	# counting using sliding window
	seq = sequence.upper() #convert to uppercase
	counts = {kmer: 0 for kmer in generate_kmers(k)} 
	n_windows = len(seq) - k + 1 # number of windows
	# if the sequence is shorter than k, there are no windows
	if n_windows <= 0:
		return counts, 0

	for i in range(n_windows): # sliding window
		kmer = seq[i:i+k] # extract kmer
		# Only consider kmers made of standard alphabet; otherwise skip
		if all(ch in 'ACGT' for ch in kmer): # check if all characters are in the standard alphabet
			counts[kmer] = counts.get(kmer, 0) + 1 # increment count

	return counts, n_windows # return counts and number of windows

# Exercise nr. 3
def print_percentages(sequence):
	print(f"Sequence (length {len(sequence)}): {sequence}\n")

	for k in (2, 3): # for dinucleotides and trinucleotides
		counts, total = count_kmers(sequence, k)
		print(f"{k}-nucleotide (k={k}) totals: {total} windows")

		# Sort kmers alphabetically for deterministic output
		for kmer in sorted(counts.keys()):
			c = counts[kmer]
			pct = (c / total * 100) if total > 0 else 0.0
			print(f"{kmer}: {c}/{total} = {pct:.4f}%")

		print('\n')


if __name__ == '__main__':
	print_percentages(S)
