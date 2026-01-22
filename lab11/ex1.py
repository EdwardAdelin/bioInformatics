import math
import pandas as pd

# --- Input Data ---
# S1: The CpG+ (Island) Sequence
S1 = "ATCGATTCGATATCATACACGTAT"

# S2: The CpG- (Non-Island) Sequence
S2 = "CTCGACTAGTATGAAGTCCACGCTTG"

# S: The Test Sequence
S_test = "CAGGTTGGAAACGTAA"

# Nucleotides
bases = ['A', 'C', 'G', 'T']

# --- Helper Functions ---

def calculate_transition_matrix(sequence, bases):
    """
    Counts transitions and converts them to probabilities (with Laplace smoothing).
    """
    # Initialize counts with 1 (Laplace smoothing) to avoid log(0) errors
    counts = {b1: {b2: 1 for b2 in bases} for b1 in bases}
    
    # Count transitions in the sequence
    for i in range(len(sequence) - 1):
        current_n = sequence[i]
        next_n = sequence[i+1]
        counts[current_n][next_n] += 1
        
    # Convert to probabilities (divide by row totals)
    probs = {b1: {} for b1 in bases}
    for b1 in bases:
        total_transitions = sum(counts[b1].values())
        for b2 in bases:
            probs[b1][b2] = counts[b1][b2] / total_transitions
            
    return probs

def calculate_log_likelihood_matrix(plus_model, minus_model, bases):
    """
    Calculates the beta matrix: log2(P+ / P-)
    """
    llr_matrix = {b1: {} for b1 in bases}
    
    for b1 in bases:
        for b2 in bases:
            p_plus = plus_model[b1][b2]
            p_minus = minus_model[b1][b2]
            
            # Calculate Log Likelihood Ratio
            # Formula: log2(Prob_Plus / Prob_Minus)
            llr = math.log2(p_plus / p_minus)
            llr_matrix[b1][b2] = round(llr, 3)
            
    return llr_matrix

def score_sequence(sequence, llr_matrix):
    """
    Scores a new sequence by summing the LLR of its transitions.
    """
    score = 0
    print(f"\nScoring Sequence: {sequence}")
    print(f"{'Transition':<12} | {'Score'}")
    print("-" * 25)
    
    for i in range(len(sequence) - 1):
        current_n = sequence[i]
        next_n = sequence[i+1]
        step_score = llr_matrix[current_n][next_n]
        score += step_score
        print(f"{current_n} -> {next_n} : {step_score:>8}")
        
    return score

# --- Step 1 & 2: Train Models ---
# Calculate probabilities for Island (+) and Non-Island (-)
prob_plus = calculate_transition_matrix(S1, bases)
prob_minus = calculate_transition_matrix(S2, bases)

# --- Step 3: Create Log-Likelihood Matrix ---
beta_matrix = calculate_log_likelihood_matrix(prob_plus, prob_minus, bases)

# Display the Matrix nicely
print("### Calculated Log-Likelihood Matrix (beta) ###")
df_matrix = pd.DataFrame(beta_matrix).T # Transpose to have 'From' as rows
print(df_matrix)

# --- Final Step: Test the New Sequence ---
final_score = score_sequence(S_test, beta_matrix)

print("-" * 25)
print(f"Total Score: {final_score:.4f}")

if final_score > 0:
    print("Result: The sequence S belongs to a CpG ISLAND.")
else:
    print("Result: The sequence S belongs to a NON-ISLAND region.")