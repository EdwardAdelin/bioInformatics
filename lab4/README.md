Codon comparison: SARS-CoV-2 vs Influenza

Files produced/used
- `NC_045512.2.fasta` — SARS‑CoV‑2 genome (already present)
- `influenza_raw.fasta` — concatenate influenza segments (you must download)
- `codon_compare.py` — analysis script
- `requirements.txt` — minimal Python deps

How to get an influenza FASTA (example using PowerShell + curl):

```
# Download segments for Influenza A (example: A/Puerto_Rico/8/1934) from NCBI
# You can download individual segment FASTAs and concatenate them into influenza_raw.fasta
curl -o influenza_seg1.fasta "https://www.ncbi.nlm.nih.gov/search/api/sequence/NC_002016.1/?report=fasta"
curl -o influenza_seg2.fasta "https://www.ncbi.nlm.nih.gov/search/api/sequence/NC_002017.1/?report=fasta"
# ... repeat for other segments or use the accession you prefer
Get-Content influenza_seg*.fasta | Set-Content influenza_raw.fasta
```

Then install deps and run:

```
python -m pip install -r requirements.txt
python codon_compare.py
```

Outputs: three PNG files in this folder:
- `covid_top10_codons.png`
- `flu_top10_codons.png`
- `top_codons_comparison.png`
