"""Use AI to: Adapt your current algorithm in order to make an app that takes a FASTA file 
and reads the seq. content from it and display the relative procentages for 
the symbols present in the alphabet of seq. 
note: FASTA represents a file format that contains DNA, ARN or proteins seq. . 
Thus it contains the information for your input."""

# Create an example FASTA file
fasta_content = """>>sequence B
ggtaagtgctctagtacaaacacccccaatattgtgatataattaaaattatattcatat
tctgttgccagattttacacttttaggctatattagagccatcttctttgaagcgttgtc
tatgcatcgatcgacgactg
"""
with open("example.fasta", "w") as f:
    f.write(fasta_content)

# Function to read a FASTA file and return the sequence
def read_fasta(filename):
    sequence = ""
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                continue  # Skip header lines
            sequence += line
    return sequence

# Read sequence from the FASTA file
seq = read_fasta("example.fasta")

# Calculate relative percentages
alf = [seq[0]]
counts = {seq[0]: 1}

for i in range(1, len(seq)):
    if seq[i] in alf:
        counts[seq[i]] += 1
    else:
        alf.append(seq[i])
        counts[seq[i]] = 1

print("Alphabet found in sequence:", alf)
print("Absolute counts:", counts)

total = sum(counts.values())
print("Relative frequencies (%):")
for key in counts:
    counts[key] = counts[key] / total * 100
print(counts)
