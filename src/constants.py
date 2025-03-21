import numpy as np

class DNAConstants:

  NUCLEOTIDES = np.fromiter("ACGT", dtype="U1")

  UNSIMILARITY_MATRIX = np.array([
        [0, 1, 1, 1],
        [1, 0, 1, 1],
        [1, 1, 0, 1],
        [1, 1, 1, 0],
    ], dtype=int)


