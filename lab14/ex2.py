import urllib.request
import matplotlib.pyplot as plt
import time

# --- 1. NATIVE DATA DOWNLOADER ---
def download_genome(accession_id, filename):
    """
    Downloads a genome from NCBI in FASTA format using native urllib.
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    url = f"{base_url}?db=nuccore&id={accession_id}&rettype=fasta&retmode=text"
    
    print(f"Downloading {accession_id} from NCBI...")
    try:
        with urllib.request.urlopen(url) as response:
            data = response.read().decode("utf-8")
            with open(filename, "w") as f:
                f.write(data)
        print(f"Saved to {filename}")
        return True
    except Exception as e:
        print(f"Error downloading {accession_id}: {e}")
        return False

# --- 2. NATIVE FASTA PARSER ---
def load_fasta(filename):
    """
    Reads a FASTA file and returns the sequence as a single string.
    Removes header lines and newlines.
    """
    with open(filename, "r") as f:
        lines = f.readlines()
    
    # Skip the first line (header) and join the rest
    # Filter out newlines to get a pure sequence string
    sequence = "".join([line.strip() for line in lines if not line.startswith(">")])
    return sequence.upper()

# --- 3. THE "IN-BETWEEN LAYER": WINDOWED COMPARISON ---
def windowed_alignment_scan(seq1, seq2, window_size=20, step=10, threshold=0.6):
    """
    This is the 'In-Between Layer' solution.
    Instead of aligning 30k x 30k bp cell-by-cell (which crashes native Python),
    we step through 'big regions' (windows) and compare them.
    
    Args:
        seq1, seq2: The DNA strings.
        window_size: Size of the chunk to compare.
        step: How much to shift the window (stride).
        threshold: Similarity percentage required to record a match (0.0 to 1.0).
        
    Returns:
        x_coords, y_coords: Arrays of coordinates where similarity > threshold.
    """
    x_coords = []
    y_coords = []
    
    len1 = len(seq1)
    len2 = len(seq2)
    
    print(f"Comparing genomes... (Seq1: {len1}bp, Seq2: {len2}bp)")
    print("This may take a moment...")
    
    # Iterate over Seq1 (Influenza)
    for i in range(0, len1 - window_size, step):
        chunk1 = seq1[i : i + window_size]
        
        # Iterate over Seq2 (COVID-19)
        for j in range(0, len2 - window_size, step):
            chunk2 = seq2[j : j + window_size]
            
            # --- Native Score Calculation (Hamming-like) ---
            matches = 0
            for k in range(window_size):
                if chunk1[k] == chunk2[k]:
                    matches += 1
            
            similarity = matches / window_size
            
            # If the region is similar enough, store the coordinate
            if similarity >= threshold:
                x_coords.append(i) # Influenza position
                y_coords.append(j) # Covid position

    return x_coords, y_coords

# --- 4. VISUALIZATION ---
def plot_alignment(x, y, label1, label2):
    plt.figure(figsize=(10, 10))
    plt.scatter(x, y, s=1, c='blue', alpha=0.5)
    
    plt.title(f"Genome Similarity Visualization\n(Windowed Local Alignment)", fontsize=14)
    plt.xlabel(f"{label1} Position (bp)", fontsize=12)
    plt.ylabel(f"{label2} Position (bp)", fontsize=12)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    
    # Force square aspect ratio for correct geometric interpretation
    plt.axis('equal')
    plt.tight_layout()
    plt.show()

# --- MAIN EXECUTION ---
if __name__ == "__main__":
    # NCBI Accession IDs
    # Influenza A (H1N1) Segment 1 (Polymerase PB2) - Representative
    INFLUENZA_ID = "NC_002023.1" 
    # SARS-CoV-2 (Complete Genome)
    COVID_ID = "NC_045512.2"     
    
    # 1. Download Data
    if download_genome(INFLUENZA_ID, "influenza.fasta") and \
       download_genome(COVID_ID, "covid19.fasta"):
        
        # 2. Load Sequences
        seq_flu = load_fasta("influenza.fasta")
        seq_cov = load_fasta("covid19.fasta")
        
        # 3. Perform the "In-Between" Layer Alignment
        # Note: We use a small window (15bp) to catch short motif similarities
        # between these very different viruses.
        start_time = time.time()
        
        # Tuning parameters for visibility:
        # Window=15, Step=10, Threshold=0.7 (70% match in that window)
        # This acts as a heuristic filter to find local alignments.
        x_hits, y_hits = windowed_alignment_scan(
            seq_flu, 
            seq_cov, 
            window_size=15, 
            step=15,    # Increasing step speeds up native code significantly
            threshold=0.65 
        )
        
        print(f"Analysis complete in {time.time() - start_time:.2f} seconds.")
        print(f"Found {len(x_hits)} regions of local similarity.")
        
        # 4. Visualize
        if len(x_hits) > 0:
            plot_alignment(x_hits, y_hits, "Influenza A (Seg 1)", "SARS-CoV-2")
        else:
            print("No significant similarity found with current settings.")