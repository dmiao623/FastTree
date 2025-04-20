import newick
from alignment import Alignment
import random

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

def random_joining(alignment: Alignment) -> newick.Node:
    """
    Construct a tree out of the given alignment by randomly joining nodes.

    Args:
        alignment (Alignment): a multiple alignment containing the sequences

    Returns:
        newick.Node: a tree, the output of the random joining
    """

    # For simplicity, we'll replace all the labels with ints.
    profiles = alignment.profile_dict
    labels = list(profiles.keys())
    N = len(labels)
    assert N > 0
    if N == 1:
        return newick.Node(labels[0])

    nodes = [newick.Node(i) for i in labels]
    active_ids = set(range(N))
    while len(active_ids) > 2:
        i, j = random.sample(list(active_ids), 2)
        # join nodes i and j
        new_node = newick.Node(descendants=[nodes[i], nodes[j]])
        active_ids.remove(i)
        active_ids.remove(j)
        active_ids.add(len(nodes))
        nodes.append(new_node)

    # join the remaining two nodes
    i, j = active_ids
    nodes[j].add_descendant(nodes[i])

    logging.info("Random-joining algorithm completed")
    return nodes[j]
