import numpy as np

from typing import Optional, Tuple, Union

from profile import Profile, profile_distance_uncorrected, profile_weighted_join
from sequence import Sequence, sequence_distance_uncorrected

class NodeInfo:
    """Stores sequence and metadata that encapsulates either a raw sequence (as a leaf node) or an
    aggregated profile (as an internal node). 

    Attributes:

        _sequence (Optional[Sequence]): The Sequence object if the node represents a leaf.

        _profile (Optional[Profile]): The Profile object if the node represents an internal node.

        _up_distance (float): The distance from this node upward to its parent node.

        _variance (float): The variance associated with this node; typically computed during tree
            construction.

        _label (Optional[str]): The label associated with this node.

    Args:

        sequence_or_profile (Union[Sequence, Profile]): The input value, which determines whether 
            the node is a sequence or a profile.

        up_distance (float, optional): The upward distance value for this node. Defaults to 0.0.

        variance (float, optional): The initial variance associated with this node. Defaults to 0.0.
    """

    def __init__(
            self,
            sequence_or_profile: Union[Sequence, Profile],
            up_distance: float = 0.0,
            variance: float = 0.0,
            label: Optional[str] = None,
        ):
        
        if isinstance(sequence_or_profile, Sequence):
            self._sequence = sequence_or_profile
            self._profile = None
        else:
            self._sequence = None
            self._profile = sequence_or_profile

        self._up_distance = up_distance
        self._variance = variance

        self._label = label

    @property
    def sequence(self) -> Optional[Sequence]: return self._sequence

    @property
    def profile(self) -> Optional[Profile]: return self._profile

    @property
    def up_distance(self) -> float: return self._up_distance

    @property
    def variance(self) -> float: return self._variance

    @property
    def label(self) -> Optional[str]: return self._label

    def set_variance(self, variance: float): self._variance = variance

def nodeinfo_distance(n1: NodeInfo, n2: NodeInfo) -> float:
    """Computes the distance between two NodeInfos.

    A node can be either a Sequence (leaf) or a Profile (internal node).

    Args:
        n1 (NodeInfo): The first node to compare.
        n2 (NodeInfo): The second node to compare.

    Returns:
        float: The uncorrected distance between the two nodes.
    """

    if n1.sequence is not None and n2.sequence is not None:
        delta = sequence_distance_uncorrected(n1.sequence, n2.sequence)
    else:
        p1 = (n1.profile if n1.profile is not None 
                else Profile.from_aligned_sequence(n1.sequence))
        p2 = (n2.profile if n2.profile is not None 
                else Profile.from_aligned_sequence(n2.sequence))
        delta = profile_distance_uncorrected(p1, p2)

    return delta - n1.up_distance - n2.up_distance

def nodeinfo_join(n1: NodeInfo, n2: NodeInfo, d: Optional[float] = None) -> Tuple[NodeInfo, float, float]:
    """Joins two NodeInfos into a single NodeInfo.

    A node can be either a Sequence (leaf) or a Profile (internal node). Returned
    NodeInfo will have a profile that represents the combined profiles / sequences
    of n1 and n2 and have correctly computed up_distance and variance.

    Args:
        n1 (NodeInfo): The first node to compare.
        n2 (NodeInfo): The second node to compare.

    Returns:
        A Tuple containing a NodeInfo object with parameters specified above,
            and two floating point integers representing the distance to the
            left and right children.
    """


    if d is None:
        d = nodeinfo_distance(n1, n2)
    v1 = n1.variance
    v2 = n2.variance

    alpha = np.clip(0.5 + (v2 - v1) / (2 * (v1 + v2)), 0, 1)
    left_dist = alpha * d
    right_dist = (1.-alpha) * d

    up_distance = (d / 2.) + abs(v1 - v2) / (2. * d)
    variance = alpha ** 2 * v1 + (1.-alpha)**2 * v2

    p1 = (n1.profile if n1.profile is not None 
            else Profile.from_aligned_sequence(n1.sequence))
    p2 = (n2.profile if n2.profile is not None 
            else Profile.from_aligned_sequence(n2.sequence))   

    p = profile_weighted_join(p1, p2, alpha, 1.-alpha)
    return (NodeInfo(p, up_distance, variance), left_dist, right_dist)
