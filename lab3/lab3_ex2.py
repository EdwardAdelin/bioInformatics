#EX2
"""
Use the AI to design an app that uses the sliding window 
methodology in order to scan a DNA sequence from a 
fasta file, and display the melting temperature 
along the sequence, by using a chart.
The chart shall have 2 signals, one for each formula.
Note: the sliding window should have 9 positions.
Make a GUI for the app.
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

		self.canvas = None

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
			self.plot_tm(positions, tm1, tm2)
		except Exception as e:
			messagebox.showerror("Error", str(e))

	def plot_tm(self, positions, tm1, tm2):
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

if __name__ == "__main__":
	root = tk.Tk()
	app = TmApp(root)
	root.mainloop()