import numpy as np

class DNAConstants:

    MONOMERS = np.fromiter("ACGT", dtype="U1")

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

"""
[TODO]: Loads one of DNAConstants or ProteinConstants
"""

UNSIMILARITY_MATRIX = DNAConstants.UNSIMILARITY_MATRIX

CORRECTION = DNAConstants.CORRECTION
