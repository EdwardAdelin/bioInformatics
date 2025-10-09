# Find in Seq S only the dinucleotides and trinucleotides that exist, 
# without using Brute Force Engine.
# In order to achive the results, one must check this 
# combinations starting from the beginning of the sequence until the end of the sequence.

# Example: S = "ABAA"

# Function to find dinucleotides and trinucleotides in a sequence
def find_nucleotides(sequence):
    dinucleotides = set()
    trinucleotides = set()

    # Iterate through the sequence to find dinucleotides
    for i in range(len(sequence) - 1):
        dinucleotides.add(sequence[i:i+2])

    # Iterate through the sequence to find trinucleotides
    for i in range(len(sequence) - 2):
        trinucleotides.add(sequence[i:i+3])

    return dinucleotides, trinucleotides

# Input sequence
S = "ABAA"

# Find dinucleotides and trinucleotides
dinucleotides, trinucleotides = find_nucleotides(S)

# Print results
print("Dinucleotides:", dinucleotides)
print("Trinucleotides:", trinucleotides)




