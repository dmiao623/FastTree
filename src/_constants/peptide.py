import numpy as np
import utils
import blosum as bl

ALPHABET = np.fromiter("ARNDCQEGHILKMFPSTWYV", dtype="U1")

bl45 = bl.BLOSUM(45)

def calc_unsim_matrix() -> list[list[float]]:
    unsimilarity_matrix = []
    for key_i in bl45.keys():
        if key_i in ALPHABET:
            matrix_row = []
            min_similarity = 0.0
            max_similarity = 0.0
            for key_j in bl45.keys():
                if key_j in ALPHABET:
                    if bl45[key_i][key_j] < min_similarity:
                        min_similarity = bl45[key_i][key_j]
                    if bl45[key_i][key_j] > max_similarity:
                        max_similarity = bl45[key_i][key_j]
            for key_j in bl45.keys():
                if key_j in ALPHABET:
                    if key_i == key_j:
                        matrix_row.append(0.0)
                    else:
                        matrix_row.append((bl45[key_i][key_j] - max_similarity) / (min_similarity - max_similarity))
            unsimilarity_matrix.append(matrix_row)
    return unsimilarity_matrix

UNSIMILARITY_MATRIX = np.array(calc_unsim_matrix(), dtype = float)

# https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation
CHARACTER_VECTORS = {
    'A' : utils.normalize(np.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'R' : utils.normalize(np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'N' : utils.normalize(np.array([0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'D' : utils.normalize(np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'C' : utils.normalize(np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'Q' : utils.normalize(np.array([0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'E' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'G' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'H' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'I' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'L' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'K' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0])),
    'M' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0])),
    'F' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0])),
    'P' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0])),
    'S' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0])),
    'T' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0])),
    'W' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0])),
    'Y' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0])),
    'V' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])),
    'B' : utils.normalize(np.array([0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'J' : utils.normalize(np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'Z' : utils.normalize(np.array([0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])),
    'X' : utils.normalize(np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])),
    # Note: it doesn't matter what the vectors of the gap positions are
    # because when we do distance calculations, we never weight them
    '-' : utils.normalize(np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])),
    '.' : utils.normalize(np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])),
    '*' : utils.normalize(np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]))
}

def CORRECTION(raw_dist: float) -> float:
    if raw_dist >= 1:
        return float("inf")
    return -1.3 * np.log(1 - raw_dist)

def IS_GAP(c : str) -> bool:
    return c in ['-', '.', '*']
