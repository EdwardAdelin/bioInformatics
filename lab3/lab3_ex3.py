#EX3
"""
Show the minimum & maximum values over the two signals. Also, allow the user
to set a threshold (like a filter) that is able to take into consideration
only the values above the threshold.
These values above the threshold should be shown to the user on a second chart,
as horizontal bars. Thus, the chunks of the signal that are above the threshold
are shown as a horizontal bar over the sequence.
Wherever the signal is below the threshold, the chart should show empty space.
"""

import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def read_fasta(filepath):
	with open(filepath, 'r') as f:
		lines = f.readlines()
	seq = ''
	for line in lines:
		if not line.startswith('>'):
			seq += line.strip().upper()
	return seq

def calc_tm_formula1(window):
	a = window.count('A')
	t = window.count('T')
	g = window.count('G')
	c = window.count('C')
	return 2 * (a + t) + 4 * (g + c)

def calc_tm_formula2(window):
	a = window.count('A')
	t = window.count('T')
	g = window.count('G')
	c = window.count('C')
	total = a + t + g + c
	if total == 0:
		return 0
	return 64.9 + 41 * (g + c - 16.4) / total

def sliding_window_tm(seq, window_size=9):
	tm1 = []
	tm2 = []
	positions = []
	for i in range(len(seq) - window_size + 1):
		window = seq[i:i+window_size]
		tm1.append(calc_tm_formula1(window))
		tm2.append(calc_tm_formula2(window))
		positions.append(i + window_size//2)
	return positions, tm1, tm2

class TmApp:
	def __init__(self, master):
		self.master = master
		master.title("DNA Melting Temperature - Sliding Window")

		self.frame = tk.Frame(master)
		self.frame.pack(padx=10, pady=10)

		self.load_btn = tk.Button(self.frame, text="Load FASTA File", command=self.load_fasta)
		self.load_btn.pack()

		# Threshold input
		self.threshold_label = tk.Label(self.frame, text="Threshold:")
		self.threshold_label.pack()
		self.threshold_entry = tk.Entry(self.frame)
		self.threshold_entry.pack()
		self.threshold_entry.insert(0, "0")

		# Min/max display
		self.stats_label = tk.Label(self.frame, text="")
		self.stats_label.pack()

		self.canvas = None
		self.canvas2 = None

	def load_fasta(self):
		filepath = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa"), ("All files", "*.*")])
		if not filepath:
			return
		try:
			seq = read_fasta(filepath)
			if len(seq) < 9:
				messagebox.showerror("Error", "Sequence too short for sliding window.")
				return
			positions, tm1, tm2 = sliding_window_tm(seq, window_size=9)

			# Show min/max values
			min1, max1 = min(tm1), max(tm1)
			min2, max2 = min(tm2), max(tm2)
			stats = f"Formula 1: min={min1}, max={max1} | Formula 2: min={min2:.2f}, max={max2:.2f}"
			self.stats_label.config(text=stats)

			# Get threshold
			try:
				threshold = float(self.threshold_entry.get())
			except ValueError:
				threshold = 0

			self.plot_tm(positions, tm1, tm2, threshold)
		except Exception as e:
			messagebox.showerror("Error", str(e))

	def plot_tm(self, positions, tm1, tm2, threshold):
		fig, ax = plt.subplots(figsize=(7,4))
		ax.plot(positions, tm1, label="Formula 1: 2(A+T)+4(G+C)")
		ax.plot(positions, tm2, label="Formula 2: 64.9+41(G+C-16.4)/N")
		ax.set_xlabel("Position (center of window)")
		ax.set_ylabel("Melting Temperature (°C)")
		ax.set_title("Melting Temperature Along DNA Sequence")
		ax.legend()
		if self.canvas:
			self.canvas.get_tk_widget().pack_forget()
		self.canvas = FigureCanvasTkAgg(fig, master=self.frame)
		self.canvas.draw()
		self.canvas.get_tk_widget().pack()

		# Second chart: horizontal bars for values above threshold
		fig2, ax2 = plt.subplots(figsize=(7,2))
		# Choose which signal to filter (here, Formula 1)
		bar_values = [v if v > threshold else None for v in tm1]
		# Plot horizontal bars
		for i, v in enumerate(bar_values):
			if v is not None:
				ax2.barh(0, 1, left=positions[i]-0.5, height=0.5, color='orange')
		ax2.set_xlim(min(positions), max(positions))
		ax2.set_yticks([])
		ax2.set_xlabel("Position (center of window)")
		ax2.set_title(f"Signal Chunks Above Threshold ({threshold})")
		if self.canvas2:
			self.canvas2.get_tk_widget().pack_forget()
		self.canvas2 = FigureCanvasTkAgg(fig2, master=self.frame)
		self.canvas2.draw()
		self.canvas2.get_tk_widget().pack()

if __name__ == "__main__":
	root = tk.Tk()
	app = TmApp(root)
	root.mainloop()