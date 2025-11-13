import random
from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt


Entrez.email = "example@example.com"
# Fetch from NCBI 
handle = Entrez.efetch(db="nucleotide", id="NC_045512.2", rettype="fasta", retmode="text")
record = SeqIO.read(handle, "fasta")
full_seq = str(record.seq)

# Random substring 
seq_length = random.randint(1000, 3000)
start_pos = random.randint(0, len(full_seq) - seq_length)
seq = full_seq[start_pos:start_pos + seq_length]
print(f"Selected sequence length: {len(seq)} nucleotides")

# 10 random samples from this sequence
samples = []
for _ in range(10):
    sample_start = random.randint(0, len(seq) - 100)
    sample_length = random.randint(100, min(3000, len(seq) - sample_start))
    sample = seq[sample_start:sample_start + sample_length]
    samples.append(sample)
print(f"Number of samples: {len(samples)}")
for i, sample in enumerate(samples):
    print(f"Sample {i+1}: length {len(sample)} nucleotides")

# Simulate electrophoresis migration based on length
lengths = [len(s) for s in samples]
sorted_indices = sorted(range(len(lengths)), key=lambda i: lengths[i])  # sort by length ascending

# Visual representation using matplotlib
fig, ax = plt.subplots(figsize=(8, 6))
y_positions = range(1, len(samples) + 1)

for i, idx in enumerate(sorted_indices):
    length = lengths[idx]
    # Draw band: thicker for visibility
    ax.plot([0, 1], [y_positions[i], y_positions[i]], 'k-', linewidth=8, solid_capstyle='butt')
    # Label with length
    ax.text(1.05, y_positions[i], f'{length} nt', va='center')

ax.set_xlim(0, 1.5)
ax.set_ylim(0, len(samples) + 1)
ax.set_title('Electrophoresis Gel Simulation\n(Short fragments migrate farther)')
ax.set_xlabel('Migration Distance (arbitrary units)')
ax.set_ylabel('Bands')
ax.set_yticks([])  # Remove y ticks
plt.tight_layout()
plt.show()
