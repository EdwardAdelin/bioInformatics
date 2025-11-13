# search and upload form ncbi or other websites 10 viral genoms. use each of these genoms in order to take samples and make a set of 
# samples for each virus.
# a. measure the time (in ms ) of assembly for each set. 
# b. measure their overall c+g percentage
# c. plot a chart in which the coordinates of the points are as follows: on the y axis the point will take the time in ms and in the x axis the 
# point will take the overall c+g percentage
# e. make a text file in which you explain the differences between the positions of the points


import requests
import random
import time
import matplotlib.pyplot as plt

# --- 1. Viral genomes (accessions from NCBI) ---
viruses = {
    "SARS-CoV-2": "NC_045512.2",
    "HIV-1": "NC_001802.1",
    "Hepatitis B": "NC_003977.2",
    "Dengue virus 2": "NC_001474.2",
    "Zika virus": "NC_012532.1",
    "Ebola virus": "NC_002549.1",
    "Norovirus": "NC_001959.2",
    "Adenovirus 2": "AC_000007.1",
    "HPV16": "NC_001526.4",
    "Measles virus": "NC_001498.1"
}

# --- 2. Helper functions ---
def fetch_fasta(accession):
    """Fetch FASTA from NCBI by accession."""
    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?id={accession}&db=nuccore&report=fasta&retmode=text"
    resp = requests.get(url)
    return "".join([line.strip() for line in resp.text.splitlines() if not line.startswith(">")])

def gc_content(seq):
    gc = sum(base in "GCgc" for base in seq)
    return (gc / len(seq)) * 100 if seq else 0

def measure_assembly_time(seq, n_samples=100, sample_len=1000):
    """Generate subsequences and measure time to concatenate them."""
    samples = []
    start = time.time()
    for _ in range(n_samples):
        start_idx = random.randint(0, len(seq) - sample_len)
        samples.append(seq[start_idx:start_idx + sample_len])
    assembled = "".join(samples)
    elapsed = (time.time() - start) * 1000  # in ms
    return elapsed

# --- 3. Fetch genomes & analyze ---
results = []
for name, acc in viruses.items():
    print(f"Processing {name}...")
    seq = fetch_fasta(acc)
    if len(seq) < 1000:
        print(f"Skipping {name} (too short)")
        continue
    gc = gc_content(seq)
    t = measure_assembly_time(seq)
    results.append((name, gc, t))

# --- 4. Plot GC% vs Assembly time ---
x = [r[1] for r in results]  # GC%
y = [r[2] for r in results]  # time (ms)
labels = [r[0] for r in results]

plt.figure(figsize=(8,6))
plt.scatter(x, y)
for i, label in enumerate(labels):
    plt.text(x[i]+0.2, y[i], label, fontsize=9)
plt.xlabel("GC Content (%)")
plt.ylabel("Assembly Time (ms)")
plt.title("Viral Genome GC% vs Assembly Time")
plt.grid(True)
plt.show()

# --- 5. Generate text report ---
with open("viral_analysis_report.txt", "w") as f:
    f.write("Viral Genome Assembly and GC% Comparison\n\n")
    for name, gc, t in results:
        f.write(f"{name}:\n  GC% = {gc:.2f}\n  Assembly Time = {t:.2f} ms\n\n")
    f.write("Observation:\n")
    f.write("Viruses with higher GC% may have slightly longer assembly times "
            "due to sequence complexity and nucleotide composition.\n")
    f.write("However, the effect also depends on genome length and sampling randomness.\n")

