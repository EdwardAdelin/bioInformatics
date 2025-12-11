from Bio import Entrez, SeqIO
from Bio.Restriction import EcoRI
import matplotlib.pyplot as plt


Entrez.email = "example@example.com" 
# Search for 10 influenza A virus complete genomes
handle = Entrez.esearch(db="nucleotide", term="influenza A virus complete genome NOT segment", retmax=10)
record = Entrez.read(handle)
ids = record["IdList"]

print(f"Found {len(ids)} genome IDs: {ids}")

# Fetch sequences
sequences = []
for id in ids:
    try:
        handle = Entrez.efetch(db="nucleotide", id=id, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(handle, "fasta")
        sequences.append(seq_record)
        print(f"Fetched sequence for {id}: {seq_record.description}")
    except Exception as e:
        print(f"Error fetching {id}: {e}")

print(f"Successfully fetched {len(sequences)} sequences")

# Digest each sequence with EcoRI
fragments_list = []
for seq_record in sequences:
    fragments = EcoRI.catalyse(seq_record.seq)
    fragment_lengths = [len(frag) for frag in fragments]
    fragments_list.append(fragment_lengths)

# Count number of fragments for each genome
num_fragments = [len(frags) for frags in fragments_list]
print(f"Number of fragments per genome: {num_fragments}")

# Find the genome with the most fragments
max_fragments = max(num_fragments)
max_idx = num_fragments.index(max_fragments)
print(f"Genome {max_idx + 1} (ID: {ids[max_idx]}) has the most DNA sequences (fragments): {max_fragments}")

# Plot the electrophoresis gel simulations
fig, axes = plt.subplots(2, 5, figsize=(15, 10))
axes = axes.flatten()

for i, frags in enumerate(fragments_list):
    # Sort fragments in descending order (larger fragments migrate slower)
    frags_sorted = sorted(frags, reverse=True)
    # Plot as horizontal bars to simulate gel bands
    y_pos = range(len(frags_sorted))
    axes[i].barh(y_pos, frags_sorted, height=0.8, color='blue', alpha=0.7)
    axes[i].set_title(f'Genome {i+1}\n{num_fragments[i]} fragments')
    axes[i].set_xlabel('Fragment Length (bp)')
    axes[i].set_ylabel('Band Position')

plt.tight_layout()
plt.savefig('electrophoresis_gel_simulation.png')
plt.show()
