import math
import re

# --- 1. The Corpus ---
eminescu_text = "cobori în jos luceafăr blând alunecând pe-o rază pătrunde-n casă și în gând și viața-mi luminează"
stanescu_text = "leoaică tânără iubirea mi-ai sărit în faţă mă pândise-n încordare mai demult"

# The suspect text (Mihai)
mihai_text = "cobori în jos luceafăr blând pe un server vechi eu scriu coduri în noapte și beau cafea rece leoaică tânără iubirea mi-ai sărit în faţă dar eu am instalat un antivirus puternic"

# --- 2. Build Transition Models ---
def build_transition_counts(text):
    words = text.split()
    counts = {}
    vocab = set(words)
    
    # Initialize simple counts
    for i in range(len(words) - 1):
        curr_w, next_w = words[i], words[i+1]
        if curr_w not in counts: counts[curr_w] = {}
        if next_w not in counts[curr_w]: counts[curr_w][next_w] = 0
        counts[curr_w][next_w] += 1
    return counts, vocab

def get_probability(w1, w2, model_counts, vocab_size):
    # Laplace Smoothing (Add-1)
    count_w1_w2 = model_counts.get(w1, {}).get(w2, 0) + 1
    count_w1_total = sum(model_counts.get(w1, {}).values()) + vocab_size
    return count_w1_w2 / count_w1_total

# Train models
em_counts, em_vocab = build_transition_counts(eminescu_text)
st_counts, st_vocab = build_transition_counts(stanescu_text)
total_vocab_size = len(em_vocab.union(st_vocab))

# --- 3. The Sliding Window Scan ---
mihai_words = mihai_text.split()
window_size = 3 # Look at a small phrase at a time

print(f"{'Text Segment (Window)':<30} | {'Score':<8} | {'Verdict'}")
print("-" * 60)

results = []

for i in range(len(mihai_words) - 1):
    w1, w2 = mihai_words[i], mihai_words[i+1]
    
    # Calculate Probabilities
    p_eminescu = get_probability(w1, w2, em_counts, total_vocab_size)
    p_stanescu = get_probability(w1, w2, st_counts, total_vocab_size)
    
    # Log Likelihood Ratio
    # If p_eminescu == p_stanescu (both unknown), log(1) = 0
    llr = math.log2(p_eminescu / p_stanescu)
    
    # Determine Verdict
    if llr > 0.5: verdict = "EMINESCU"
    elif llr < -0.5: verdict = "STANESCU"
    else: verdict = "ORIGINAL (Mihai)"
    
    print(f"{w1 + ' ' + w2:<30} | {llr:6.2f}   | {verdict}")