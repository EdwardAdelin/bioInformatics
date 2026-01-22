# 1.use a random dna sequence of about 50 letters. use this sequence to compute the transition probabilities between letters. 
# your output should be the transition matrix stored as a json file. 
# 2.use a random english text of about 300 letters (that implies spaces, punctuation etc) and compile the transition probabilities 
# between words. store the transition matrix as a json file. for ease of implementation, you can represent each new word by using 
# a symbol of your choice (ascii).

import random
import json
import string

def compute_transition_matrix(sequence, available_states=None):
    transitions = {}
    
    for i in range(len(sequence) - 1):
        current_item = sequence[i]
        next_item = sequence[i+1]
        
        if current_item not in transitions:
            transitions[current_item] = {}
        
        if next_item not in transitions[current_item]:
            transitions[current_item][next_item] = 0
            
        transitions[current_item][next_item] += 1

    all_states = available_states if available_states else list(transitions.keys())
    
    probability_matrix = {}
    
    for state in all_states:
        probability_matrix[state] = {}
        
        if state in transitions:
            total_transitions = sum(transitions[state].values())
            for next_state, count in transitions[state].items():
                probability_matrix[state][next_state] = round(count / total_transitions, 4)
        else:
            probability_matrix[state] = {}

    return probability_matrix

def run_dna_task():
    print("--- Generating DNA Sequence ---")
    bases = ["A", "C", "G", "T"]
    dna_seq = "".join(random.choices(bases, k=50))
    print(f"Sequence: {dna_seq}")
    
    matrix = compute_transition_matrix(dna_seq, available_states=bases)
    
    filename = 'dna_transitions.json'
    with open(filename, 'w') as f:
        json.dump(matrix, f, indent=4)
    print(f"Saved: {filename}\n")

def run_text_task():
    print("--- Generating English Text ---")
    text_corpus = (
        "The sun dipped below the horizon, painting the sky in shades of orange and pink. "
        "The wind whispered through the trees, carrying the scent of pine and damp earth. "
        "A small bird chirped in the distance, calling out to its mate before the night settled in. "
        "Stars began to twinkle, one by one, like tiny diamonds scattered across a vast, dark velvet canvas. "
        "It was a peaceful evening, a moment of quiet reflection before the world woke up again."
    )
    
    translator = str.maketrans('', '', string.punctuation)
    clean_text = text_corpus.translate(translator).lower()
    words = clean_text.split()
    
    print(f"Text Length: {len(text_corpus)} chars")
    print(f"Word Count: {len(words)} words")

    unique_words = sorted(list(set(words)))
    
    word_to_symbol = {word: chr(33 + i) for i, word in enumerate(unique_words)}
    symbol_to_word = {v: k for k, v in word_to_symbol.items()}
    
    symbol_sequence = [word_to_symbol[w] for w in words]
    
    matrix = compute_transition_matrix(symbol_sequence, available_states=list(word_to_symbol.values()))
    
    output_data = {
        "legend_mapping": symbol_to_word,
        "transition_matrix": matrix
    }
    
    filename = 'word_transitions.json'
    with open(filename, 'w') as f:
        json.dump(output_data, f, indent=4)
    print(f"Saved: {filename}")

if __name__ == "__main__":
    run_dna_task()
    run_text_task()