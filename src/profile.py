import numpy as np

from numpy.typing import NDArray
from typing import Union

import constants
from alignment import Alignment
from sequence import Sequence
from constants import ALPHALEN

class Profile:
    """A class representing a profile matrix.

    [ToDo]: Update documentation and class clarification.

    Properties:

        profile (NDArray[NDArray[float]]): the profile matrix

        profile_length (int): the length of the profile

        num_sequences (int): the number of sequences stored

        ungapped (NDArray[float]): the proportion of non-gaps in each column
    """

    def __init__(self, alignment_or_sequence: Union[Alignment, Sequence, str]):
        if isinstance(alignment_or_sequence, Alignment):
            alignment = alignment_or_sequence

            s_len = alignment.alignment_length
            profile = np.zeros((ALPHALEN, s_len), dtype=float)
            for i in range(ALPHALEN):
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

        elif isinstance(alignment_or_sequence, Sequence):
            sequence = alignment_or_sequence

            s_len = sequence.sequence_length
            profile = np.zeros((ALPHALEN, s_len), dtype=float)
            
            for col_idx, char_idx in enumerate(sequence):
                if char_idx != -1:
                    profile[char_idx, col_idx] = 1.0

            self._profile = profile
            self._profile_length = s_len

        elif isinstance(alignment_or_sequence, str):
            # edgar dijkstra rolling in his grave rn

            sequence = alignment_or_sequence
            s_len = len(sequence)
        
            profile = np.zeros((constants.ALPHALEN, s_len), dtype=float)
            ungapped = np.zeros(s_len, dtype=float)
            
            for col_idx, char in enumerate(sequence):
                try:
                    profile[:, col_idx] = constants.NUCLEIC_ACID_VECTORS[char]
                    if not constants.IS_GAP(char):
                        ungapped[col_idx] = 1
                except KeyError:
                    raise ValueError(f"Encountered unknown character: {char}")

            self._profile = profile
            self._profile_length = s_len
            self._num_sequences = 1
            self._ungapped = ungapped

        # [ToDo]: Fix Me!
        self._num_sequences = 0
        self._ungapped = np.zeros(s_len, dtype=float)


    @property
    def profile(self) -> NDArray[float]: return self._profile

    @property
    def profile_length(self) -> int: return self._profile_length

    @property
    def num_sequences(self) -> int: return self._num_sequences

    @property
    def ungapped(self) -> NDArray[float]: return self._ungapped


def profile_distance_uncorrected(p1: Profile, p2: Profile) -> float:
    """Compute the distance between two profiles.

    Specifically, profile distance is calculated by
      1. Summing the position-wise dissimilarity
      2. Weighting by the fraction of non-gapped position in each profile
    
    [ToDo]: add eigendecomposition optimization

    Args:

        p1 (Profile): The first profile to compare.

        p2 (Profile): The second profile to compare.

    Returns:
        float: The corrected distance between the two profiles.
               - If all columns are gapped (no overlap), returns 0.0.
               - If the raw distance is at or beyond a model limit, returns float("inf").
    """

    p1_mat, p2_mat = p1.profile, p2.profile
    profile_length = p1.profile_length
    
    pos_dissimilarities = np.array(
        [p1_mat[:, i] @ (constants.UNSIMILARITY_MATRIX @ p2_mat[:, i])
        for i in range(profile_length)],
    dtype=float)
    
    column_weights = p1.ungapped * p2.ungapped

    column_mask = column_weights > 0.0
    if not np.any(column_mask):
        return 0.0

    raw_dist = np.sum(
        pos_dissimilarities[column_mask] * column_weights[column_mask]) \
        / np.sum(column_weights[column_mask])

    return raw_dist


def profile_distance_corrected(p1: Profile, p2: Profile) -> float:
    """Computes the corrected distance between two profiles.

    Specifically, profile distance is calculated by
      1. Summing the position-wise dissimilarity
      2. Weighting by the fraction of non-gapped position in each profile
      3. Applies a correction to compensate for hidden substitutions

    Args:

        p1 (Profile): The first profile to compare.

        p2 (Profile): The second profile to compare.

    Returns:
        float: The corrected distance between the two profiles.
               - If all columns are gapped (no overlap), returns 0.0.
               - If the raw distance is at or beyond a model limit, returns float("inf").
    """

    raw_dist = profile_distance_uncorrected(p1, p2)
    corrected_dist = constants.CORRECTION(raw_dist)
    return corrected_dist
