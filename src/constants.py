import numpy as np

class DNAConstants:

    ALPHABET = np.fromiter("ACGT", dtype="U1")

    UNSIMILARITY_MATRIX = np.array([
            [0, 1, 1, 1],
            [1, 0, 1, 1],
            [1, 1, 0, 1],
            [1, 1, 1, 0],
        ], dtype=int)

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

CORRECTION = DNAConstants.CORRECTION
