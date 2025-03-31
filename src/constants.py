import numpy as np

class DNAConstants:

    MONOMERS = np.fromiter("ACGT", dtype="U1")

    UNSIMILARITY_MATRIX = np.array([
            [0, 1, 1, 1],
            [1, 0, 1, 1],
            [1, 1, 0, 1],
            [1, 1, 1, 0],
        ], dtype=int)

    CORRECTION = lambda d : -0.75 * np.log(1 - (4/3)*d)
    CORRECTION_LIMIT = 0.75

"""
[TODO]: Loads one of DNAConstants or ProteinConstants
"""

UNSIMILARITY_MATRIX = DNAConstants.UNSIMILARITY_MATRIX

CORRECTION = 
