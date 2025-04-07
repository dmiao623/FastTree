import numpy as np

class DNAConstants:

    ALPHABET = np.fromiter("ACGT", dtype="U1")

    UNSIMILARITY_MATRIX = np.array([
            [0, 1, 1, 1],
            [1, 0, 1, 1],
            [1, 1, 0, 1],
            [1, 1, 1, 0],
        ], dtype=int)

    def _normalize(A: np.typing.NDArray) -> np.typing.NDArray:
        s = np.sum(A)
        if s == 0:
            return A
        return A / s

    # https://en.wikipedia.org/wiki/FASTA_format#Sequence_representation
    NUCLEIC_ACID_VECTORS = {
        'A' : _normalize(np.array([1, 0, 0, 0])),
        'C' : _normalize(np.array([0, 1, 0, 0])),
        'G' : _normalize(np.array([0, 0, 1, 0])),
        'T' : _normalize(np.array([0, 0, 0, 1])),
        'U' : _normalize(np.array([0, 0, 0, 1])),
        'R' : _normalize(np.array([1, 0, 1, 0])),
        'Y' : _normalize(np.array([0, 1, 0, 1])),
        'K' : _normalize(np.array([0, 0, 1, 1])),
        'M' : _normalize(np.array([1, 1, 0, 0])),
        'S' : _normalize(np.array([0, 1, 1, 0])),
        'W' : _normalize(np.array([1, 0, 0, 1])),
        'B' : _normalize(np.array([0, 1, 1, 1])),
        'D' : _normalize(np.array([1, 0, 1, 1])),
        'H' : _normalize(np.array([1, 1, 0, 1])),
        'V' : _normalize(np.array([1, 1, 1, 0])),
        'N' : _normalize(np.array([1, 1, 1, 1])),
        '-' : _normalize(np.array([0, 0, 0, 0])),
        '.' : _normalize(np.array([0, 0, 0, 0])),
    }

    def CORRECTION(raw_dist: float) -> float:
        if raw_dist >= 0.75:
            return float("inf")
        return -0.75 * np.log(1 - (4/3)*raw_dist)

class ProteinConstants:
    """
    [ToDo]: add protein constants
    """

    pass

"""
[TODO]: Loads one of DNAConstants or ProteinConstants
"""

ALPHABET = DNAConstants.ALPHABET

UNSIMILARITY_MATRIX = DNAConstants.UNSIMILARITY_MATRIX

NUCLEIC_ACID_VECTORS = DNAConstants.NUCLEIC_ACID_VECTORS

CORRECTION = DNAConstants.CORRECTION
