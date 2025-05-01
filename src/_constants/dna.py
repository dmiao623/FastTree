import numpy as np
import utils

ALPHABET = np.fromiter("ACGT", dtype="U1")

UNSIMILARITY_MATRIX = np.array([
        [0, 1, 1, 1],
        [1, 0, 1, 1],
        [1, 1, 0, 1],
        [1, 1, 1, 0],
    ], dtype=int)

# https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation
CHARACTER_VECTORS = {
    'A' : utils.normalize(np.array([1, 0, 0, 0])),
    'C' : utils.normalize(np.array([0, 1, 0, 0])),
    'G' : utils.normalize(np.array([0, 0, 1, 0])),
    'T' : utils.normalize(np.array([0, 0, 0, 1])),
    'U' : utils.normalize(np.array([0, 0, 0, 1])),
    'R' : utils.normalize(np.array([1, 0, 1, 0])),
    'Y' : utils.normalize(np.array([0, 1, 0, 1])),
    'K' : utils.normalize(np.array([0, 0, 1, 1])),
    'M' : utils.normalize(np.array([1, 1, 0, 0])),
    'S' : utils.normalize(np.array([0, 1, 1, 0])),
    'W' : utils.normalize(np.array([1, 0, 0, 1])),
    'B' : utils.normalize(np.array([0, 1, 1, 1])),
    'D' : utils.normalize(np.array([1, 0, 1, 1])),
    'H' : utils.normalize(np.array([1, 1, 0, 1])),
    'V' : utils.normalize(np.array([1, 1, 1, 0])),
    'N' : utils.normalize(np.array([1, 1, 1, 1])),
    # Note: it doesn't matter what the vectors of the gap positions are
    # because when we do distance calculations, we never weight them
    '-' : utils.normalize(np.array([1, 1, 1, 1])),
    '.' : utils.normalize(np.array([1, 1, 1, 1])),
}

def CORRECTION(raw_dist: float) -> float:
    if raw_dist >= 0.75:
        return float("inf")
    return -0.75 * np.log(1 - (4/3)*raw_dist)

def IS_GAP(c : str) -> bool:
    return c in ['-', '.']
