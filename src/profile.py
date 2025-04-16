import numpy as np

from numpy.typing import NDArray
from typing import Union

import constants
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

    def __init__(self, p_mat: NDArray[NDArray[float]], num_sequences: int, ungapped: NDArray[float]):
        self._profile = p_mat
        self._profile_length = len(p_mat[0])
        self._num_sequences = num_sequences
        self._ungapped = ungapped

    @classmethod
    def from_aligned_sequence(self, aligned_seq: str):
        s_len = len(aligned_seq)
    
        profile = np.zeros((constants.ALPHALEN, s_len), dtype=float)
        ungapped = np.zeros(s_len, dtype=float)
        
        for col_idx, char in enumerate(aligned_seq):
            try:
                profile[:, col_idx] = constants.NUCLEIC_ACID_VECTORS[char]
                if not constants.IS_GAP(char):
                    ungapped[col_idx] = 1
            except KeyError:
                raise ValueError(f"Encountered unknown character: {char}")

        return Profile(profile, 1, ungapped)

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

def profile_weighted_join(p1: Profile, p2: Profile, w1: float, w2: float) -> Profile:
    """
    Compute the profile formed by joining profiles p1 and p2 with weights w1 and w2.
    """
    if w1 == 0. and w2 == 0.:
        w1, w2 = 1., 1.
    freq_mat_1 = p1.profile * p1.ungapped * p1.num_sequences
    freq_mat_2 = p2.profile * p2.ungapped * p2.num_sequences
    freq_mat_weighted = freq_mat_1 * w1 + freq_mat_2 * w2
    counts = np.sum(freq_mat_weighted, axis=0)
    new_p_mat = np.divide(freq_mat_weighted,
                          counts,
                          out=(0.25 * np.ones_like(freq_mat_weighted)),
                          where=(counts != 0))
    new_ungapped = counts / (p1.num_sequences * w1 + p2.num_sequences * w2)

    return Profile(new_p_mat,
                   p1.num_sequences + p2.num_sequences,
                   new_ungapped)
