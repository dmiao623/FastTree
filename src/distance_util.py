import numpy as np

from typing import TypeVar

import constants
from profile import Profile

Sequence = TypeVar("Sequence", bound=AbstractSequence)
Node = Union[Sequence, Profile]

def sequence_distance(s1: Sequence, s2: Sequence) -> float:
    """Computes the distance between two biological sequences

    Args:
        s1 (AbstractSequence): the first sequence to compare
        s2 (AbstractSequence): the second sequence to compare

    Returns:
        float: The distance between the two sequences.

    Raises:
        ValueError: Raised if provided sequences have different lengths.
    """

    if s1.sequence_length != s2.sequence_length:
        raise ValueError("Sequences must be of the same length to compute distance.")
    mismatches = np.sum(seq1.sequence != seq2.sequence)

    raw_dist = mismatches / seq1.sequence.size

    corrected_dist = constants.CORRECTION(raw_dist)

    return corrected_dist

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
        float: The corrected distance between the two profiles.
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

    corrected_dist = constants.CORRECTION(raw_dist)

    return corrected_dist

def node_distance(n1: Node, n2: Node) -> float:
    """Computes the distance between two nodes.

    A node can be either a Sequence (leaf) or a Profile (internal node).

    Args:
        n1 (Node): the first node to compare
        n2 (Node): the second node to compare

    Returns:
        float: The correct distance between the two nodes.
    """

    if isinstance(n1, Sequence) and isinstance(n2, Sequence):
        return sequence_distance(n1, n2)
    
    p1 = n1 if isinstance(n1, Profile) else Profile.from_sequence(n1)
    p2 = n2 if isinstance(n2, Profile) else Profile.from_sequence(n2)
    return profile_distance(p1, p2)
