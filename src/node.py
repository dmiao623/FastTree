import numpy as np

from numpy.typing import NDArray
from typing import Optional, Tuple, Union

from profile import Profile, profile_distance_uncorrected
from sequence import Sequence, sequence_distance_uncorrected

class Node:
    def __init__(
            self,
            sequence_or_profile: Union[Sequence, Profile],
            up_distance: float = 0.0,
            variance: float = 0.0,
            left_child: Optional[Tuple[float, "Node"]] = None,
            right_child: Optional[Tuple[float, "Node"]] = None,
        ):
        
        if isinstance(sequence_or_profile, Sequence):
            self._sequence = sequence_or_profile
            self._profile = None
        else:
            self._sequence = None
            self._profile = sequence_or_profile

        self._up_distance = up_distance
        self._variance = variance

        if left_child is not None:
            self._left_child_dist = left_child[0]
            self._left_child = left_child[1]
        else: 
            self._left_child_dist = None
            self._left_child = None

        if right_child is not None:
            self._right_child_dist = right_child[0]
            self._right_child = right_child[1]
        else: 
            self._right_child_dist = None
            self._right_child = None

    @property
    def sequence(self) -> Optional[Sequence]: return self._sequence

    @property
    def profile(self) -> Optional[Profile]: return self._profile

    @property
    def up_distance(self) -> float: return self._up_distance

    @property
    def variance(self) -> float: return self._variance

def node_distance(n1: Node, n2: Node) -> float:
    """Computes the distance between two nodes.

    A node can be either a Sequence (leaf) or a Profile (internal node).

    Args:
        n1 (Node): the first node to compare
        n2 (Node): the second node to compare

    Returns:
        float: The correct distance between the two nodes.
    """

    if n1.sequence is not None and n2.sequence is not None:
        delta = sequence_distance_uncorrected(n1.sequence, n2.sequence)
    else:
        p1 = (n1.profile if n1.profile is not None 
                else Profile.from_sequence(n1.sequence))
        p2 = (n2.profile if n2.profile is not None 
                else Profile.from_sequence(n2.sequence))
        delta = profile_distance_uncorrected(p1, p2)

    return delta - n1.up_distance - n2.up_distance

def node_join(n1: Node, n2: Node, d: Optional[float] = None) -> Node:
    if d is None:
        d = node_distance(n1, n2)
    v1 = n1.variance
    v2 = n2.variance

    up_distance = (d / 2.) + abs(v1 - v2) / (2. * d)
    variance = (v1 + v2) / 4. + ((v1 - v2) / (2. * d)) ** 2

    left_dist = 0.5 * d + (abs(v1 - v2) / (2. * d))
    right_dist = d - left_dist

    left_dist = max(0., left_dist)
    right_dist = max(0., right_dist)

    w_left = left_dist / d
    w_right = right_dist / d
    p1 = (n1.profile if n1.profile is not None 
            else Profile.from_sequence(n1.sequence))
    p2 = (n2.profile if n2.profile is not None 
            else Profile.from_sequence(n2.sequence))   

    p_mat = p1.profile * w_left + p2.profile * w_right
    p = Profile(p_mat)

    return Node(p, up_distance, variance, (left_dist, v1), (right_dist, v2))
