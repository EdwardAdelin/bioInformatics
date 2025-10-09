"""FASTA sliding-window frequency analyzer with GUI.

This script provides a small Tkinter application that lets the user
choose a FASTA file from disk, select a sliding-window size (default 30),
and plots the relative frequencies of symbols (e.g., A,C,G,T for DNA)
computed on every sliding window across the sequence.

Usage: run this file with a Python interpreter that has tkinter and matplotlib
installed (both are commonly available). The application embeds a Matplotlib
figure inside a Tkinter window.

Features:
- Open FASTA files (single-sequence FASTA supported)
- Choose sliding window size
- Show sequence name and length
- Plot relative frequencies along the sequence (one line per symbol)

This file is intended as a self-contained example for bioinformatics labs.
"""
"""Design an application using AI, which 
contains a GUI which allows the user to choose
 a FASTA file. The content of the file should be analysed by 
 using a sliding window of 30 positions. The content for 
 each sliding window should be used in order to extract 
 the relative freqs of the symbols found in the alphabet 
 of the seq. Thus, your input should be the DNA sequence 
 from the fasta file and the output should be the values 
 of the relative freqs of each symbol from the seq 
 translated as lines on a chart. Thus, your chart in the 
 case of DNA should have 4 lines which reflect the values 
 found over the seq"""
from __future__ import annotations

import os
from typing import Dict, List, Tuple

# Attempt to import matplotlib with TkAgg backend for GUI embedding. If that
# fails (for example Tk isn't available on the system), fall back to the
# non-interactive Agg backend and run in a simple CLI/headless mode that
# writes a PNG file instead of showing a GUI.
MATPLOTLIB_HEADLESS = False
try:
	import matplotlib
	matplotlib.use("TkAgg")
	from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
	import matplotlib.pyplot as plt
except Exception:
	# try non-interactive backend
	try:
		import matplotlib
		matplotlib.use("Agg")
		import matplotlib.pyplot as plt
		MATPLOTLIB_HEADLESS = True
	except Exception as e:  # pragma: no cover - import failure
		raise ImportError("matplotlib is required to run this app (either with TkAgg for GUI or Agg for headless plotting)") from e

if not MATPLOTLIB_HEADLESS:
	import tkinter as tk
	from tkinter import filedialog, messagebox


def parse_fasta(path: str) -> Tuple[str, str]:
	"""Parse a simple FASTA file and return (header, sequence).

	Only the first sequence in the file is returned. Lines starting with
	'>' are considered headers. Sequence lines are concatenated and returned
	in uppercase with all whitespace removed.
	"""
	if not os.path.exists(path):
		raise FileNotFoundError(path)
	header = ""
	seq_lines: List[str] = []
	with open(path, "r", encoding="utf-8") as fh:
		for line in fh:
			line = line.strip()
			if not line:
				continue
			if line.startswith(">"):
				if header:
					# already read one sequence; stop (single-sequence mode)
					break
				header = line[1:].strip()
			else:
				seq_lines.append(line)
	sequence = "".join(seq_lines).upper()
	return header or os.path.basename(path), sequence


def sliding_window_freqs(sequence: str, window: int) -> Tuple[List[int], Dict[str, List[float]]]:
	"""Compute relative frequencies per sliding window.

	Returns (positions, freqs) where positions are the center indices (1-based)
	of each window and freqs is a dict mapping symbol->list of relative
	frequencies (floats between 0 and 1) for each window.
	"""
	seq = sequence.replace("\n", "").replace(" \t", "")
	n = len(seq)
	if window <= 0:
		raise ValueError("window must be > 0")
	if n < window:
		return [], {}

	# Determine alphabet from sequence in a stable order: A,C,G,T then others
	base_order = ["A", "C", "G", "T"]
	unique = sorted(set(seq))
	alphabet = [b for b in base_order if b in unique] + [u for u in unique if u not in base_order]

	counts: Dict[str, List[float]] = {sym: [] for sym in alphabet}
	positions: List[int] = []

	for i in range(0, n - window + 1):
		window_seq = seq[i : i + window]
		positions.append(i + window // 2 + 1)  # 1-based center position
		total = len(window_seq)
		# count each symbol
		freq_map: Dict[str, int] = {}
		for ch in window_seq:
			freq_map[ch] = freq_map.get(ch, 0) + 1
		for sym in alphabet:
			counts[sym].append(freq_map.get(sym, 0) / total)

	return positions, counts


class FastaFreqApp:
	def __init__(self, master: tk.Tk) -> None:
		self.master = master
		master.title("FASTA sliding-window frequency analyzer")
		master.geometry("900x600")

		# Top controls
		frm = tk.Frame(master)
		frm.pack(side=tk.TOP, fill=tk.X, padx=8, pady=8)

		self.open_btn = tk.Button(frm, text="Open FASTA...", command=self.open_fasta)
		self.open_btn.pack(side=tk.LEFT)

		tk.Label(frm, text="Window size:").pack(side=tk.LEFT, padx=(10, 2))
		self.win_var = tk.StringVar(value="30")
		self.win_entry = tk.Entry(frm, width=6, textvariable=self.win_var)
		self.win_entry.pack(side=tk.LEFT)

		self.analyze_btn = tk.Button(frm, text="Analyze & Plot", command=self.analyze_and_plot)
		self.analyze_btn.pack(side=tk.LEFT, padx=(10, 0))

		self.info_label = tk.Label(master, text="No file loaded", anchor="w")
		self.info_label.pack(fill=tk.X, padx=8)

		# Matplotlib figure
		self.fig, self.ax = plt.subplots(figsize=(8, 5))
		self.canvas = FigureCanvasTkAgg(self.fig, master=master)
		self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=1)

		self.filepath: str | None = None
		self.sequence: str = ""
		self.header: str = ""

	def open_fasta(self) -> None:
		path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fa *.fasta *.fna"), ("All files", "*")])
		if not path:
			return
		try:
			header, seq = parse_fasta(path)
		except Exception as e:
			messagebox.showerror("Error", f"Could not read FASTA: {e}")
			return
		self.filepath = path
		self.header = header
		self.sequence = seq
		self.info_label.config(text=f"Loaded: {os.path.basename(path)} — {len(seq)} nt — header: {header}")
		# clear previous plot
		self.ax.clear()
		self.ax.set_title("(click 'Analyze & Plot' to compute)")
		self.canvas.draw()

	def analyze_and_plot(self) -> None:
		if not self.sequence:
			messagebox.showinfo("No file", "Please open a FASTA file first.")
			return
		try:
			window = int(self.win_var.get())
		except ValueError:
			messagebox.showerror("Invalid window", "Window size must be an integer.")
			return

		positions, freqs = sliding_window_freqs(self.sequence, window)
		if not positions:
			messagebox.showwarning("Too short", "Sequence is shorter than the sliding window.")
			return

		self.ax.clear()
		# colors: use matplotlib default cycle
		for sym, values in freqs.items():
			self.ax.plot(positions, values, label=sym)

		self.ax.set_xlabel("Position (center of window)")
		self.ax.set_ylabel("Relative frequency")
		self.ax.set_ylim(-0.02, 1.02)
		self.ax.set_title(f"Symbol relative frequencies — window={window} — {os.path.basename(self.filepath or '')}")
		self.ax.legend(title="Symbols", bbox_to_anchor=(1.05, 1), loc="upper left")
		self.ax.grid(True, linestyle="--", alpha=0.4)
		self.fig.tight_layout()
		self.canvas.draw()


def main() -> None:
	if MATPLOTLIB_HEADLESS:
		# Headless mode: simple CLI flow. Ask user for a FASTA file via input
		print("Running in headless mode (no Tk). Provide a FASTA file path.")
		path = input("FASTA path (leave empty to exit): ").strip()
		if not path:
			print("No file provided. Exiting.")
			return
		try:
			header, seq = parse_fasta(path)
		except Exception as e:
			print(f"Error reading FASTA: {e}")
			return
		try:
			window = int(input("Window size [30]: ") or "30")
		except ValueError:
			print("Invalid window size. Using 30.")
			window = 30
		positions, freqs = sliding_window_freqs(seq, window)
		if not positions:
			print("Sequence shorter than window; nothing to plot.")
			return

		fig, ax = plt.subplots(figsize=(10, 4))
		for sym, values in freqs.items():
			ax.plot(positions, values, label=sym)
		ax.set_xlabel("Position (center of window)")
		ax.set_ylabel("Relative frequency")
		ax.set_ylim(-0.02, 1.02)
		ax.set_title(f"Symbol relative frequencies — window={window} — {os.path.basename(path)}")
		ax.legend(title="Symbols")
		ax.grid(True, linestyle="--", alpha=0.4)
		fig.tight_layout()
		out = os.path.splitext(os.path.basename(path))[0] + f"_window{window}.png"
		fig.savefig(out)
		print(f"Plot saved to {out}")
	else:
		root = tk.Tk()
		app = FastaFreqApp(root)
		root.mainloop()


if __name__ == "__main__":
	main()

