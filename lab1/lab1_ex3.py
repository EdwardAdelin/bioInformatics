"""Use AI to: Adapt your current algorithm in order to make an app that takes a FASTA file 
and reads the seq. content from it and display the relative procentages for 
the symbols present in the alphabet of seq. 
note: FASTA represents a file format that contains DNA, ARN or proteins seq. . 
Thus it contains the information for your input."""

# Create an example FASTA file in the same folder as this script
import os
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EXAMPLE_FASTA = os.path.join(SCRIPT_DIR, "example.fasta")

fasta_content = ">>sequence B\nggtaagtgctctagtacaaacacccccaatattgtgatataattaaaattatattcatat\ntctgttgccagattttacacttttaggctatattagagccatcttctttgaagcgttgtc\ntatgcatcgatcgacgactg\n"
with open(EXAMPLE_FASTA, "w") as f:
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

# Refactor: implement the original counting algorithm in its own function
def compute_frequencies(seq):
    """Return (alphabet_list, absolute_counts_dict, relative_percentages_dict)
    Uses the same algorithm/ordering as the original script (first-seen ordering).
    """
    if not seq:
        return [], {}, {}

    alf = [seq[0]]
    counts = {seq[0]: 1}

    for i in range(1, len(seq)):
        if seq[i] in alf:
            counts[seq[i]] += 1
        else:
            alf.append(seq[i])
            counts[seq[i]] = 1

    total = sum(counts.values())
    percentages = {}
    for key in counts:
        percentages[key] = counts[key] / total * 100

    return alf, counts, percentages

# GUI using tkinter
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

class FastaApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FASTA Alphabet Frequencies")
        self.geometry("600x400")

        # Default to the example FASTA located next to this script
        self.filename_var = tk.StringVar(value=EXAMPLE_FASTA)

        top_frame = tk.Frame(self)
        top_frame.pack(fill=tk.X, padx=8, pady=8)

        tk.Label(top_frame, text="FASTA file:").pack(side=tk.LEFT)
        tk.Entry(top_frame, textvariable=self.filename_var, width=40).pack(side=tk.LEFT, padx=6)
        tk.Button(top_frame, text="Browse...", command=self.browse_file).pack(side=tk.LEFT)
        tk.Button(top_frame, text="Load & Analyze", command=self.load_and_analyze).pack(side=tk.LEFT, padx=6)

        # Treeview for results
        cols = ("symbol", "absolute", "relative")
        self.tree = ttk.Treeview(self, columns=cols, show="headings", height=10)
        self.tree.heading("symbol", text="Symbol")
        self.tree.heading("absolute", text="Absolute Count")
        self.tree.heading("relative", text="Relative (%)")
        self.tree.column("symbol", width=100, anchor=tk.CENTER)
        self.tree.column("absolute", width=120, anchor=tk.CENTER)
        self.tree.column("relative", width=120, anchor=tk.CENTER)
        self.tree.pack(fill=tk.BOTH, expand=True, padx=8, pady=(0,8))

        bottom_frame = tk.Frame(self)
        bottom_frame.pack(fill=tk.X, padx=8, pady=4)
        self.alphabet_label = tk.Label(bottom_frame, text="Alphabet: ")
        self.alphabet_label.pack(side=tk.LEFT)

        # Load default example file on start
        self.after(100, lambda: self.load_and_analyze())

    def browse_file(self):
        fname = filedialog.askopenfilename(title="Open FASTA file", filetypes=[("FASTA files", "*.fasta;*.fa;*.txt"), ("All files", "*.*")])
        if fname:
            self.filename_var.set(fname)

    def load_and_analyze(self):
        fname = self.filename_var.get()
        if not fname:
            messagebox.showwarning("No file", "Please choose a FASTA file.")
            return
        if not os.path.exists(fname):
            messagebox.showerror("File not found", f"File does not exist: {fname}")
            return

        seq = read_fasta(fname)
        if not seq:
            messagebox.showinfo("Empty sequence", "No sequence data found in the selected FASTA file.")
            return

        alf, counts, percentages = compute_frequencies(seq)

        # Update alphabet label
        self.alphabet_label.config(text=f"Alphabet: {', '.join(alf)}")

        # Update treeview
        for row in self.tree.get_children():
            self.tree.delete(row)

        for symbol in alf:
            abs_count = counts[symbol]
            rel = percentages[symbol]
            self.tree.insert("", tk.END, values=(symbol, abs_count, f"{rel:.2f}"))

        # Also print to console for compatibility with original script behavior
        print("Alphabet found in sequence:", alf)
        print("Absolute counts:", counts)
        print("Relative frequencies (%):")
        print({k: f"{v:.2f}" for k, v in percentages.items()})


if __name__ == "__main__":
    app = FastaApp()
    app.mainloop()
