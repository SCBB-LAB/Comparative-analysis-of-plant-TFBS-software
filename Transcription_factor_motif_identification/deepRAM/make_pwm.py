from Bio import SeqIO
import numpy as np
import sys
def sequences_to_pwm(sequences):
    sequence_length = len(sequences[0])
    pwm_matrix = np.zeros((4, sequence_length))
    for seq in sequences:
        for i, base in enumerate(seq):
            if base == "A":
                pwm_matrix[0, i] += 1
            elif base == "C":
                pwm_matrix[1, i] += 1
            elif base == "G":
                pwm_matrix[2, i] += 1
            elif base == "T":
                pwm_matrix[3, i] += 1
    # Normalize the PWM
    pwm_matrix = pwm_matrix / len(sequences)
    return pwm_matrix

def pwm_to_jaspar(pwm_matrix, matrix_id="M0001", matrix_name="ExampleMatrix"):
    if pwm_matrix.shape[0] != 4:
        raise ValueError("The PWM matrix must have 4 rows corresponding to A, C, T, G.")    
    jaspar_lines = ["\t".join(map(str, pwm_matrix[0])), "\t".join(map(str, pwm_matrix[1])), "\t".join(map(str, pwm_matrix[3])), "\t".join(map(str, pwm_matrix[2]))]    
    return "\n".join(jaspar_lines)

fasta_file = sys.argv[1] + "_logo.fa"
sequences = [str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")]

pwm_result = sequences_to_pwm(sequences)
jaspar_format = pwm_to_jaspar(pwm_result)
output_file = sys.argv[1] + "_pwm.txt"
with open(output_file, "w") as file:
    file.write(jaspar_format)


