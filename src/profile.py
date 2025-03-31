import numpy as np

from numpy.typing import NDArray

from alignment import Alignment
import constants

class Profile:
    def __init__(self, alignment: Alignment):
        a_len = alignment.alphabet_size
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
        self._alphabet_length = a_len
        self._profile_length = s_len

    @property
    def profile(self) -> NDArray[float]: return self._profile

    @property
    def profile_length(self) -> int: return self._profile_length

def profile_distance(p1: Profile, p2: Profile) -> float:
    """Compute the distance between two profiles.

    Specifically, profile distance is calculated by
      1. Summing the position-wise dissimilarity
      2. Weighting by the fraction of non-gapped position in each profile
      3. Applies a correction to compensate for hidden substitutions
    
    [ToDo]: add eigendecomposition optimization

    Args:
        p1 (Profile): the first profile to compare
        p2 (Profile): the second profile to compare

    Returns:
        float: The corrected distance between these two profiles.
               - If all columns are gapped (no overlap), returns 0.0.
               - If the raw distance is at or beyond a model limit, returns float('inf').

    """

    p1_mat, p2_mat = p1.profile, p2.profile
    profile_length = p1.profile_length
    
    pos_dissimilarities = np.array(
        [p1_mat[:, i] @ (constants.UNSIMILARITY_MATRIX @ p2_mat[:, i])
        for i in range(profile_length)],
    dtype=float)
    
    p1_nongap = np.sum(p1_mat, axis=0)
    p2_nongap = np.sum(p2_mat, axis=0)
    column_weights = p1_nongap * p2_nongap

    column_mask = column_weights > 0.0
    if not np.any(column_mask):
        return 0.0

    raw_dist = np.sum(
        pos_dissimilarities[column_mask] * column_weights[column_mask]) \
        / np.sum(column_weights[column_mask])

    if raw_dist >= constants.CORRECTION_LIMIT:
        return float("inf")
    else:
        corrected_dist = constants.CORRECTION(raw_dist)

    return corrected_dist

