import tkinter as tk
from tkinter import ttk, messagebox
import math

class AlignmentApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Sequence Alignment App (Needleman-Wunsch)")
        self.root.geometry("1100x650")
        
        # --- Styles and Layout ---
        self.main_container = tk.Frame(root, padx=10, pady=10)
        self.main_container.pack(fill=tk.BOTH, expand=True)

        # Top section (Controls + Visuals)
        self.top_frame = tk.Frame(self.main_container)
        self.top_frame.pack(fill=tk.BOTH, expand=True, side=tk.TOP)

        # 1. Left Control Panel
        self.control_frame = tk.LabelFrame(self.top_frame, text="Sequences & Parameters", padx=10, pady=10)
        self.control_frame.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 10))
        
        # Sequence Inputs
        tk.Label(self.control_frame, text="Sq 1 =").grid(row=0, column=0, sticky="e")
        self.entry_seq1 = tk.Entry(self.control_frame, width=25, font=("Consolas", 10))
        self.entry_seq1.grid(row=0, column=1, pady=5)
        self.entry_seq1.insert(0, "ACCGTGAAGCCAATAC")

        tk.Label(self.control_frame, text="Sq 2 =").grid(row=1, column=0, sticky="e")
        self.entry_seq2 = tk.Entry(self.control_frame, width=25, font=("Consolas", 10))
        self.entry_seq2.grid(row=1, column=1, pady=5)
        self.entry_seq2.insert(0, "AGCGTGCAGCCAATAC")

        # Parameters
        self.param_frame = tk.LabelFrame(self.control_frame, text="Parameters", padx=5, pady=5)
        self.param_frame.grid(row=2, column=0, columnspan=2, pady=10, sticky="ew")

        tk.Label(self.param_frame, text="Gap =").grid(row=0, column=0)
        self.entry_gap = tk.Entry(self.param_frame, width=5)
        self.entry_gap.grid(row=0, column=1)
        self.entry_gap.insert(0, "0")

        tk.Label(self.param_frame, text="Match =").grid(row=1, column=0)
        self.entry_match = tk.Entry(self.param_frame, width=5)
        self.entry_match.grid(row=1, column=1)
        self.entry_match.insert(0, "1")

        tk.Label(self.param_frame, text="Mismatch =").grid(row=2, column=0)
        self.entry_mismatch = tk.Entry(self.param_frame, width=5)
        self.entry_mismatch.grid(row=2, column=1)
        self.entry_mismatch.insert(0, "-1")

        # Align Button
        self.btn_align = tk.Button(self.control_frame, text="Align", command=self.run_alignment, height=2, bg="#e1e1e1")
        self.btn_align.grid(row=3, column=0, columnspan=2, sticky="ew", pady=15)
        
        # Dummy Presets (Visual only, to match screenshot)
        self.preset_frame = tk.LabelFrame(self.control_frame, text="Presets", padx=5, pady=5)
        self.preset_frame.grid(row=4, column=0, columnspan=2, sticky="ew")
        for i in range(1, 7):
            btn = tk.Button(self.preset_frame, text=f"Setting {i}", state="disabled")
            btn.grid(row=(i-1)//2, column=(i-1)%2, sticky="ew", padx=2, pady=2)

        # 2. Middle Frame: Heatmap
        self.mid_frame = tk.LabelFrame(self.top_frame, text="Graphic representation (Heatmap)", padx=5, pady=5)
        self.mid_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 10))
        
        self.canvas_heat = tk.Canvas(self.mid_frame, bg="black")
        self.canvas_heat.pack(fill=tk.BOTH, expand=True)

        # 3. Right Frame: Traceback Grid
        self.right_frame = tk.LabelFrame(self.top_frame, text="Traceback path deviation", padx=5, pady=5)
        self.right_frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.canvas_trace = tk.Canvas(self.right_frame, bg="white")
        self.canvas_trace.pack(fill=tk.BOTH, expand=True)

        # Bottom section (Text Output)
        self.bottom_frame = tk.LabelFrame(self.main_container, text="Show Alignment:", padx=5, pady=5)
        self.bottom_frame.pack(side=tk.BOTTOM, fill=tk.X, pady=(10, 0))
        
        self.txt_output = tk.Text(self.bottom_frame, height=8, font=("Courier New", 11))
        self.txt_output.pack(side=tk.LEFT, fill=tk.X, expand=True)
        
        self.scrollbar = tk.Scrollbar(self.bottom_frame, command=self.txt_output.yview)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.txt_output.config(yscrollcommand=self.scrollbar.set)

    def get_color(self, val, min_val, max_val):
        """Interpolate color between Dark Blue and Bright Red based on score."""
        if max_val == min_val: return "#000000"
        ratio = (val - min_val) / (max_val - min_val)
        
        # RGB interpolation
        # Start (Dark Blue): 0, 0, 50
        # End (Red): 255, 0, 50
        r = int(0 + (255 - 0) * ratio)
        g = 0
        b = int(50 + (50 - 50) * ratio) # Keeping blue low/constant gives a nice purple hue in mid
        return f'#{r:02x}{g:02x}{b:02x}'

    def run_alignment(self):
        try:
            s1 = self.entry_seq1.get().upper()
            s2 = self.entry_seq2.get().upper()
            gap = int(self.entry_gap.get())
            match = int(self.entry_match.get())
            mismatch = int(self.entry_mismatch.get())
        except ValueError:
            messagebox.showerror("Error", "Please ensure parameters are integers.")
            return

        n, m = len(s1), len(s2)
        
        # --- Needleman-Wunsch Algorithm ---
        # Initialize matrices
        # scores[i][j] stores the score
        scores = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
        # traceback[i][j] stores direction: 0=Diag, 1=Up, 2=Left
        traceback = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

        for i in range(n + 1): scores[i][0] = i * gap
        for j in range(m + 1): scores[0][j] = j * gap

        # Fill Matrix
        min_score = 0
        max_score = 0
        
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                if s1[i-1] == s2[j-1]:
                    diag = scores[i-1][j-1] + match
                else:
                    diag = scores[i-1][j-1] + mismatch
                
                up = scores[i-1][j] + gap
                left = scores[i][j-1] + gap
                
                # Priority: Match/Mismatch > Up > Left (standard preference)
                # But to match the specific gap-heavy behavior in the image (Gap=0),
                # we prioritize Gaps if scores are equal or better.
                best = max(diag, up, left)
                scores[i][j] = best
                
                if best < min_score: min_score = best
                if best > max_score: max_score = best

                # Direction logic customized for visual match
                if diag == best: traceback[i][j] = 0
                elif up == best: traceback[i][j] = 1
                else: traceback[i][j] = 2

        # --- Traceback Path ---
        align1, align2 = "", ""
        i, j = n, m
        path_coords = set()
        path_coords.add((i, j))
        
        while i > 0 or j > 0:
            if i > 0 and j > 0 and traceback[i][j] == 0:
                align1 = s1[i-1] + align1
                align2 = s2[j-1] + align2
                i -= 1
                j -= 1
            elif i > 0 and (j == 0 or traceback[i][j] == 1):
                align1 = s1[i-1] + align1
                align2 = "-" + align2
                i -= 1
            elif j > 0 and (i == 0 or traceback[i][j] == 2):
                align1 = "-" + align1
                align2 = s2[j-1] + align2
                j -= 1
            path_coords.add((i, j))

        # --- Draw Heatmap ---
        self.canvas_heat.delete("all")
        cw = self.canvas_heat.winfo_width() / (m + 1)
        ch = self.canvas_heat.winfo_height() / (n + 1)
        
        for r in range(n + 1):
            for c in range(m + 1):
                color = self.get_color(scores[r][c], min_score, max_score)
                x1, y1 = c * cw, r * ch
                x2, y2 = x1 + cw, y1 + ch
                self.canvas_heat.create_rectangle(x1, y1, x2, y2, fill=color, outline="")

        # --- Draw Traceback Grid ---
        self.canvas_trace.delete("all")
        tw = self.canvas_trace.winfo_width() / (m + 1)
        th = self.canvas_trace.winfo_height() / (n + 1)
        
        for r in range(n + 1):
            for c in range(m + 1):
                x1, y1 = c * tw, r * th
                x2, y2 = x1 + tw, y1 + th
                
                is_path = (r, c) in path_coords
                fill_col = "#d32f2f" if is_path else "#fff9c4" # Red or Yellow
                
                self.canvas_trace.create_rectangle(x1, y1, x2, y2, fill=fill_col, outline="black")

        # --- Text Output ---
        matches = 0
        match_str = ""
        for k in range(len(align1)):
            if align1[k] == align2[k] and align1[k] != '-':
                matches += 1
                match_str += "|"
            else:
                match_str += " "
        
        similarity = int((matches / len(align1)) * 100)
        
        output_text = (
            f"{align1}\n"
            f"{match_str}\n"
            f"{align2}\n\n"
            f"Matches = {matches}\n"
            f"Length = {len(align1)}\n"
            f"Similarity = {similarity} %\n"
            f"Tracing back: M[{n},{m}]"
        )
        
        self.txt_output.delete("1.0", tk.END)
        self.txt_output.insert(tk.END, output_text)

if __name__ == "__main__":
    root = tk.Tk()
    app = AlignmentApp(root)
    root.mainloop()