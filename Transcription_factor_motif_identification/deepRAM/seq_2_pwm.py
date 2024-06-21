from Bio import SeqIO
import numpy as np
import sys, os
def sequences_to_pwm(sequences):
    # Get the length of the motif
    motif_length = len(sequences[0])

    # Convert sequences to a list of lists (matrix)
    matrix = [list(seq) for seq in sequences]

    # Transpose the matrix to get positions as columns
    matrix_transposed = np.array(matrix).T

    # Create an empty PWM matrix
    pwm_matrix = np.zeros((4, motif_length), dtype=float)

    # Define the order of nucleotides
    nucleotides = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    # Populate the PWM matrix
    for position, column in enumerate(matrix_transposed):
        total_count = len(column)
        for nucleotide, index in nucleotides.items():
            count = np.count_nonzero(column == nucleotide)
            pwm_matrix[index, position] = count / total_count

    return pwm_matrix

def print_pwm(pwm_matrix):
    nucleotides = ['A', 'C', 'T', 'G']
    print("\t".join(nucleotides))
    for position in range(pwm_matrix.shape[1]):
        values = [f"{pwm_matrix[i, position]:.3f}" for i in range(pwm_matrix.shape[0])]
        print("\t".join(values))

if __name__ == "__main__":
    # Read sequences from a file
    input_file = sys.argv[1]
    sequences = [str(record.seq) for record in SeqIO.parse(input_file, "fasta")]

    # Convert sequences to PWM
    pwm_matrix = sequences_to_pwm(sequences)

    # Print the PWM matrix
    #print("Position Weight Matrix:")
    #print("------------------------")
    print_pwm(pwm_matrix)



