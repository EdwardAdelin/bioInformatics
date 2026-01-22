import numpy as np

def predict_n_states(transition_matrix, initial_vector, steps=5):
    """
    Predicts the state vector for a given number of discrete steps.

    Args:
        transition_matrix (np.array): The square transition matrix (NxN).
        initial_vector (np.array): The initial state vector (N).
        steps (int): Number of steps to predict.

    Returns:
        None (Prints the results directly).
    """
    # 1. Validation to ensure matrix and vector dimensions match
    rows, cols = transition_matrix.shape
    vec_len = len(initial_vector)
    
    if rows != cols:
        raise ValueError("The transition matrix must be square.")
    if cols != vec_len:
        raise ValueError(f"Matrix columns ({cols}) must match vector size ({vec_len}).")

    print(f"--- Initial State (t=0) ---\n{initial_vector}\n")

    current_vector = initial_vector.copy()

    # 2. Iterate through the steps
    for t in range(1, steps + 1):
        # Perform Matrix Multiplication: v(t+1) = M * v(t)
        # np.dot is the standard dot product for arrays
        next_vector = np.dot(transition_matrix, current_vector)
        
        print(f"--- Step {t} ---")
        print(next_vector)
        print() # Empty line for readability
        
        # Update the current vector for the next iteration
        current_vector = next_vector

# --- Example Usage: DNA Substitution Model (4-States: A, C, G, T) ---
if __name__ == "__main__":
    # Define an arbitrary 4x4 Matrix (Example: DNA mutation probabilities)
    # Rows/Cols represent A, C, G, T
    # Note: In a valid probability matrix, columns usually sum to 1.
    matrix_data = [
        [0.7, 0.1, 0.1, 0.1],  # To A
        [0.1, 0.7, 0.1, 0.1],  # To C
        [0.1, 0.1, 0.7, 0.1],  # To G
        [0.1, 0.1, 0.1, 0.7]   # To T
    ]
    
    # Define an Initial Vector (Example: 100% chance starting at 'A')
    vector_data = [1.0, 0.0, 0.0, 0.0]

    # Convert to NumPy arrays
    M = np.array(matrix_data)
    v0 = np.array(vector_data)

    print("Running Prediction Software...\n")
    predict_n_states(M, v0, steps=5)