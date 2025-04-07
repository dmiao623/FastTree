import numpy as np

from numpy.typing import NDArray

import constants
#from alignment import Alignment
from sequence import Sequence

class Profile:
    """
    A class representing a profile matrix.

    Properties:
        profile (NDArray[NDArray[float]]): the profile matrix
        profile_length (int): the length of the profile
        num_sequences (int): the number of sequences stored

    Note: because gaps are represented by the vector [0, 0, 0, 0],
    the frequency of non-gaps can be calculated using np.sum(profile, axis=0).
    """

    # TODO: fix __init__ and from_sequence methods
    #def __init__(self, alignment: Alignment):
    def __init__(self, alignment):
        s_len = alignment.alignment_length

        profile = np.zeros((a_len, s_len), dtype=float)
        for i in range(a_len):
            profile[i] = np.sum(alignment.alignment == i, axis=0)
        non_gap_counts = np.sum(alignment.alignment != -1, axis=0)
        with np.errstate(divide='ignore', invalid='ignore'):
            profile = np.where(
                non_gap_counts > 0,
                profile / non_gap_counts,
                0.0
            )

        self._profile = profile
        self._profile_length = s_len

    @classmethod
    def from_sequence(cls, sequence: Sequence) -> "Profile":
        s_len = sequence.sequence_length
        
        profile = np.zeros((constants.ALPHABET_SIZE, s_len), dtype=float)
        
        for col_idx, char_idx in enumerate(sequence):
            if char_idx != -1:
                profile[char_idx, col_idx] = 1.0

        new_profile = cls.__new__(cls)
        new_profile._profile = profile
        new_profile._profile_length = s_len

        return new_profile

    @classmethod
    def from_aligned_sequence(cls, sequence: str) -> "Profile":
        """
        Takes as input a sequence that may have gaps and non-ACGT characters
        and returns the corresponding profile.
        """
        s_len = len(sequence)
        
        profile = np.zeros((len(constants.ALPHABET), s_len), dtype=float)
        
        for col_idx, char in enumerate(sequence):
            try:
                profile[:, col_idx] = constants.NUCLEIC_ACID_VECTORS[char]
            except KeyError:
                raise ValueError(f"Encountered unknown character: {char}")

        new_profile = cls.__new__(cls)
        new_profile._profile = profile
        new_profile._profile_length = s_len
        new_profile._num_sequences = 1

        return new_profile

    @property
    def profile(self) -> NDArray[NDArray[float]]: return self._profile

    @property
    def profile_length(self) -> int: return self._profile_length

    @property
    def num_sequences(self) -> int: return self._num_sequences
