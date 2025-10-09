# Find in Seq S only the dinucleotides and trinucleotides that exist, 
# without using Brute Force Engine.
# In order to achive the results, one must check this 
# combinations starting from the beginning of the sequence until the end of the sequence.

# Example: S = "ABAA"

# Function to find dinucleotides and trinucleotides in a sequence
def find_nucleotides_with_counts(sequence):
    dinucleotides = {}
    trinucleotides = {}

    # Iterate through the sequence to find dinucleotides and their counts
    for i in range(len(sequence) - 1):
        dinucleotide = sequence[i:i+2]
        if dinucleotide in dinucleotides:
            dinucleotides[dinucleotide] += 1
        else:
            dinucleotides[dinucleotide] = 1

    # Iterate through the sequence to find trinucleotides and their counts
    for i in range(len(sequence) - 2):
        trinucleotide = sequence[i:i+3]
        if trinucleotide in trinucleotides:
            trinucleotides[trinucleotide] += 1
        else:
            trinucleotides[trinucleotide] = 1

    return dinucleotides, trinucleotides

# Input sequence
S = "ABAA"

# Find dinucleotides and trinucleotides with counts
dinucleotides, trinucleotides = find_nucleotides_with_counts(S)

# Print results
print("Dinucleotides:")
for dinucleotide, count in dinucleotides.items():
    print(f"{dinucleotide}: {count}")

print("\nTrinucleotides:")
for trinucleotide, count in trinucleotides.items():
    print(f"{trinucleotide}: {count}")




