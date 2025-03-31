import numpy as np

from numpy.typing import NDArray

import constants
from alignment import Alignment

class Profile:
    def __init__(self, alignment: Alignment):
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

    @property
    def profile(self) -> NDArray[float]: return self._profile

    @property
    def profile_length(self) -> int: return self._profile_length
