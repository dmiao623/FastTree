import newick
from profile import profile_distance_uncorrected
from alignment import Alignment
from typing import Dict

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

def neighbor_joining(alignment: Alignment) -> newick.Node:
    """
    Run the neighbor-joining algorithm, with no optimization, on the given
    alignment. (We compute the distance matrix in O(N^2L) time, then we
    calculate each of the N-2 joins in O(N^2) time, for a total time complexity
    of O(N^2(N+L).)

    Args:
        alignment (Alignment): a multiple alignment containing the sequences

    Returns:
        newick.Node: a tree, the output of the neighbor-joining algorithm
    """

    # For efficiency, we'll replace all the labels with ints.
    profiles = alignment.profile_dict
    labels = list(profiles.keys())
    N = len(labels)
    assert N > 0
    if N == 1:
        return newick.Node(labels[0])

    # First, construct the distance matrix (which will be stored as
    # a Dict[int, Dict[int, float]]).
    logger.info(f"Constructing distance matrix for {N} profiles")
    distance_matrix = dict()
    for i, l1 in enumerate(labels):
        distance_matrix[i] = dict()
        for j, l2 in enumerate(labels):
            if i == j:
                distance_matrix[i][j] = 0
            elif i < j:
                distance_matrix[i][j] = \
                    profile_distance_uncorrected(profiles[l1], profiles[l2])
            else:
                distance_matrix[i][j] = distance_matrix[j][i]
        logger.info(f"Computed distances from node {i}")
    logger.info("Successfully constructed distance matrix")

    # Next, do N-2 joins.
    logger.info("Starting neighbor-joining phase")
    edges = []
    for join in range(N-2):
        n = N - join # The current number of nodes in the tree.
        assert len(distance_matrix) == n
        assert all(len(subdict) == n for subdict in distance_matrix.values())

        # Calculate the neighbor-joining matrix.
        d_star = dict()
        total_distance = {i: sum(subdict.values())\
            for i, subdict in distance_matrix.items()}
        for i, subdict in distance_matrix.items():
            d_star[i] = dict()
            for j, dist in subdict.items():
                if i == j:
                    d_star[i][j] = 0
                else:
                    d_star[i][j] = (n - 2) * distance_matrix[i][j]\
                        - total_distance[i] - total_distance[j]
        
        # Get the minimum non-diagonal element in d_star.
        mn = (float('inf'), -1, -1)
        for i, subdict in d_star.items():
            for j, dist in subdict.items():
                if i != j:
                    mn = min(mn, (dist, i, j))
        _, mn_i, mn_j = mn

        # Calculate limb lengths.
        delta = (total_distance[mn_i] - total_distance[mn_j]) / (n - 2)
        limb_length_i = (distance_matrix[mn_i][mn_j] + delta) / 2
        limb_length_j = (distance_matrix[mn_i][mn_j] - delta) / 2

        # Add a new node m to the tree and update the distance matrix.
        # First, add m.
        m = N + join
        distance_matrix[m] = dict()
        distance_matrix[m][m] = 0
        for k in distance_matrix:
            dist = (distance_matrix[k][mn_i] + distance_matrix[k][mn_j]\
                - distance_matrix[mn_i][mn_j]) / 2
            distance_matrix[m][k] = dist
            distance_matrix[k][m] = dist
        # Then, remove mn_i and mn_j.
        del distance_matrix[mn_i]
        del distance_matrix[mn_j]
        for subdict in distance_matrix.values():
            del subdict[mn_i]
            del subdict[mn_j]

        # Add the new edges to the tree, and repeat.
        edges.append((mn_i, m, limb_length_i))
        edges.append((mn_j, m, limb_length_j))
        logger.info(f"Join {join+1} of {N-2} completed")

    # At this point, the tree has only two nodes. Join them.
    i, j = distance_matrix.keys()
    edges.append((i, j, distance_matrix[i][j]))
    root = j

    # Now, all we have to do is construct the tree in Newick format.
    nodes = [newick.Node(labels[i] if i < N else None) for i in range(2*N-2)]
    for i, j, weight in edges:
        nodes[i].length = weight
        nodes[j].add_descendant(nodes[i])

    logger.info("Neighbor-joining algorithm completed")
    return nodes[root]
