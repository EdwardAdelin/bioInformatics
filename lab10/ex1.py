import math
import pandas as pd

seqs_all = [
    "GAGGTAAAC",  # row 1
    "TCCGTAAGT",  # row 2
    "CAGGTTGGA",  # row 3
    "ACAGTCAGT",  # row 4
    "TAGGTCATT",  # row 5
    "TAGGTACTG",  # row 6
    "ATGGTAACT",  # row 7
    "CAGGTATAC",  # row 8
    "TGTGTGAGT",  # row 9
    "AAGGTAAGT",  # row 10
]

# Choose how many sequences to use (your prompt says 9, image shows 10).
USE_FIRST_N = 9   # <- change to 10 if your instructor says to use all 10
seqs = seqs_all[:USE_FIRST_N]

# Basic checks
if not seqs:
    raise ValueError("No sequences provided.")
L = len(seqs[0])
if any(len(s) != L for s in seqs):
    raise ValueError("All sequences must have the same length (aligned motifs).")
alphabet = ["A", "C", "G", "T"]
if any(any(ch not in alphabet for ch in s) for s in seqs):
    raise ValueError("Sequences must contain only A/C/G/T.")

counts = {b: [0] * L for b in alphabet}
for s in seqs:
    for j, ch in enumerate(s):
        counts[ch][j] += 1

count_df = pd.DataFrame(counts, index=range(1, L + 1)).T  # rows A,C,G,T ; cols 1..L
count_df.columns = [str(i) for i in range(1, L + 1)]

N = len(seqs)
relfreq_df = count_df / N

pseudocount = 1.0  # set to 0.0 if you want weights == relative frequencies
denom = N + pseudocount * len(alphabet)
weight_df = (count_df + pseudocount) / denom

null = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}  # uniform background

loglike_df = weight_df.copy()
for b in alphabet:
    loglike_df.loc[b] = [math.log(p / null[b]) for p in weight_df.loc[b]]

#display results
pd.set_option("display.precision", 4)

print(f"Using N={N} sequences, motif length L={L}\n")

print("COUNT MATRIX (A/C/G/T x positions)")
print(count_df, "\n")

print("RELATIVE FREQUENCIES MATRIX (counts / N)")
print(relfreq_df, "\n")

print(f"WEIGHT MATRIX (smoothed probs, pseudocount={pseudocount})")
print(weight_df, "\n")

print("LOG-LIKELIHOOD MATRIX (ln(weight / null))")
print(loglike_df, "\n")


S = "CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA"

motif_len = L  # = 9

def score_window(window, loglike_df):
    """Compute log-likelihood score of a single window"""
    score = 0.0
    for j, base in enumerate(window, start=1):
        score += float(loglike_df.loc[base, str(j)])
    return score

results = []
for i in range(len(S) - motif_len + 1):
    window = S[i:i + motif_len]
    score = score_window(window, loglike_df)
    results.append({
        "Start (1-based)": i + 1,
        "Window": window,
        "Score": score
    })

scan_df = pd.DataFrame(results)

print("SLIDING WINDOW SCORES")
print(scan_df.to_string(index=False))

print("\nHigh-scoring windows (Score > 0):")
signals = scan_df[scan_df["Score"] > 0].sort_values("Score", ascending=False)
print(signals.to_string(index=False))

print(f"\nTotal windows scanned: {len(scan_df)}")
print(f"Number of candidate signals (Score > 0): {len(signals)}")
