# use the transition matrix from the json matrix in order to synthesize new sequences of text based on the transition matrix.

import json
import random
import os

def load_matrix_from_file(filename):
    if not os.path.exists(filename):
        print(f"Error: The file '{filename}' was not found.")
        return None
        
    with open(filename, 'r') as f:
        data = json.load(f)
    return data

def synthesize_sequence(matrix, length=20, start_key=None):
    current_key = start_key
    
    if not current_key or current_key not in matrix:
        valid_starts = [k for k in matrix.keys() if matrix[k]] 
        if not valid_starts:
            return []
        current_key = random.choice(valid_starts)

    sequence = [current_key]

    for _ in range(length - 1):
        transitions = matrix.get(current_key, {})
        
        if not transitions:
            break

        next_keys = list(transitions.keys())
        probabilities = list(transitions.values())
        
        next_key = random.choices(next_keys, weights=probabilities, k=1)[0]
        
        sequence.append(next_key)
        current_key = next_key

    return sequence

def main():
    print("DNA SYNTHESIS")
    dna_data = load_matrix_from_file('dna_transitions.json')
    
    if dna_data:
        dna_seq_list = synthesize_sequence(dna_data, length=50, start_key="A")
        print(f"Generated Sequence: {''.join(dna_seq_list)}\n")

    print("TEXT SYNTHESIS")
    text_data = load_matrix_from_file('word_transitions.json')
    
    if text_data:
        word_matrix = text_data["transition_matrix"]
        legend = text_data["legend_mapping"]
        
        symbol_seq = synthesize_sequence(word_matrix, length=20, start_key="O")
        
        word_seq = [legend.get(symbol, "???") for symbol in symbol_seq]
        
        print(f"Generated Text: {' '.join(word_seq)}")

if __name__ == "__main__":
    main()