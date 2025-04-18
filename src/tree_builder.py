import math

from constants import CORRECTION
from node_info import NodeInfo, nodeinfo_distance, nodeinfo_join
from sequence import Sequence
from alignment import Alignment
from utils import UnionFind
import newick

from typing import List, Optional, Set, Tuple

NodeID = int

import logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG)

class TreeBuilder:
    """Constructs a phylogenetic tree using the FastTree algorithm.

    Attributes:

        _num_sequences (int): The number of input sequences.

        _tophits_threshold (int): The threshold factor (approximately sqrt(num_sequences))that 
            limits candidate join nodes.

        _refresh_interval (int): The interval (in steps) at which top hits are recomputed.

        _distance_cache (List[List[float]]): A triangular cache storing pairwise distances between
            nodes. (Lower triangular, so _distance_cache[i][j] is defined iff i > j)

        _nodes (List[TreeBuilder.Node]): List containing all nodes (both initial and merged) in the
            tree.

        _num_nodes (int): The total number of nodes (original sequences + merged nodes).

        _active_ids (Set[int]): Set of active node IDs that have not yet been merged.

        _steps (int): Counter for the number of join steps executed.

        _union_find (UnionFind): Union-find data structure used to efficiently manage node
            groupings during merges.
        
    Args:

        alignment (Alignment): The alignment used to construct the tree.

        thresh_cp (int, optional): A multiplier used to compute the top-hit threshold. Higher
            factor is safer but slower. Defaults to 2.

        refresh_interval (Optional[int], optional): The interval for refreshing top-hit candidates.
            If not provided, no refreshes will be performed.
    """

    class Node:
        """Internal representation of a tree node in TreeBuilder.

        This class encapsulates node information including its unique identifier, associated 
        NodeInfo, top-hit candidate nodes, and optional pointers to its left and right child nodes
        resulting from a merge.

        Attributes:

            _id (int): Unique identifier for the node.

            _node_info (NodeInfo): Associated NodeInfo that holds sequence/profile and distance
                data.

            tophit_ids (Set[int]): Set of node IDs that are considered top-hit candidates for
                potential merging.

            _leftchild_id (Optional[int]): Node ID of the left child after a node join.

            _rightchild_id (Optional[int]): Node ID of the right child after a node join.

        Args:

            id (int): Unique identifier for the node.

            node_info (NodeInfo): The NodeInfo instance containing sequence/profile and distance
                information.

            tophit_ids (Optional[Set[int]], optional): A set of candidate node IDs for joining.
                Defaults to an empty set.

            leftchild_id (Optional[int], optional): Identifier for the left child node if it exists.

            rightchild_id (Optional[int], optional): Identifier for the right child node if it
                exists.

        """

        def __init__(self,
                     id: int,
                     node_info: NodeInfo,
                     tophit_ids: Optional[Set[int]]=None,
                     leftchild_id: Optional[int]=None,
                     rightchild_id: Optional[int]=None):
            self._id = id
            self._node_info = node_info

            self.tophit_ids = tophit_ids if tophit_ids else set()

            self._leftchild_id = leftchild_id
            self._rightchild_id = rightchild_id

        @property
        def id(self) -> int: return self._id

        @property
        def node_info(self) -> NodeInfo: return self._node_info

        @property
        def leftchild_id(self) -> Optional[int]: return self._leftchild_id

        @property
        def rightchild_id(self) -> Optional[int]: return self._rightchild_id

    def __init__(self,
                 alignment: Alignment,
                 thresh_cp: int=2,
                 refresh_interval: Optional[int]=None,
                 enable_tophits_approx=True):
        logging.info("Initializing tree builder")
        self._num_sequences = alignment.alignment_size
        self._tophits_threshold = thresh_cp*math.isqrt(self._num_sequences)
        self._refresh_interval = refresh_interval if refresh_interval else 2*self._num_sequences
        self._enable_tophits_approx = enable_tophits_approx

        node_infos = [NodeInfo(profile, label=label) for label, profile in alignment.profile_dict.items()]
        self._distance_cache = [
            [-1 for j in range(i)] for i in range(self._num_sequences)
        ]
        self._nodes = [
            TreeBuilder.Node(i, node_info) for i, node_info in enumerate(node_infos)
        ]

        logging.info("Constructing top hits")
        self._num_nodes = self._num_sequences
        self._active_ids = set(range(self._num_sequences))
        self._recompute_tophits()

        self._steps = 0
        self._union_find = UnionFind(2*self._num_sequences)
        logging.info("Initialization of tree builder completed")

    def _distance_util(self, nd_id1: NodeID, nd_id2: NodeID):
        """Computes distance between two nodes identified by their IDs. Uses a cached distance
        matrix to avoid redundant computations.

        Parameters:

            nd_id1 (NodeID): Identifier of the first node.

            nd_id2 (NodeID): Identifier of the second node.

        Returns:
            float: The computed or cached distance between the two nodes.
        """
        if nd_id1 == nd_id2:
            return 0
        nd_id1, nd_id2 = max(nd_id1, nd_id2), min(nd_id1, nd_id2)
        if self._distance_cache[nd_id1][nd_id2] == -1:
            distance = nodeinfo_distance(self._nodes[nd_id1].node_info,
                                         self._nodes[nd_id2].node_info)
            self._distance_cache[nd_id1][nd_id2] = distance
        return self._distance_cache[nd_id1][nd_id2]

    def _node_join(self, nd_id1: NodeID, nd_id2: NodeID):
        """Joins two active nodes into a new parent node. The parent node is set to be active
        and the two provided nodes are set to be children.

        Parameters:

            nd_id1 (NodeID): Identifier of the first node to be merged.

            nd_id2 (NodeID): Identifier of the second node to be merged.
        """
        logging.info(f"Joining nodes {nd_id1} and {nd_id2}")

        assert nd_id1 in self._active_ids
        assert nd_id2 in self._active_ids

        nd1 = self._nodes[nd_id1]
        nd2 = self._nodes[nd_id2]

        id = self._num_nodes
        self._distance_cache.append([-1] * self._num_nodes)
        self._num_nodes += 1

        node_info = nodeinfo_join(nd1.node_info, nd2.node_info)
        self._union_find.union(id, nd_id1)
        self._union_find.union(id, nd_id2)

        potential_tophit_ids = [
            self._union_find.find(nd_id)
                for nd_id in nd1.tophit_ids | nd2.tophit_ids
                    if nd_id != id
        ]
        potential_tophit_ids = list(set(potential_tophit_ids))
        self._nodes.append(TreeBuilder.Node(id, node_info, set(), nd_id1, nd_id2))
        self._compute_single_tophits_list(id, potential_tophit_ids)

        self._active_ids.add(id)
        self._active_ids.remove(nd_id1)
        self._active_ids.remove(nd_id2)

    def _compute_single_tophits_list(self, nd_id: NodeID, candidates=None):
        """Compute the top-hits list of a single node.
        """
        logging.info(f"Computing top-hits list of node {nd_id}")
        if candidates is None:
            candidates = list(self._active_ids)
        sorted_node_ids = sorted(candidates,
                                 key=lambda j: self._distance_util(nd_id, j))
        assert sorted_node_ids[0] == nd_id
        tophits = sorted_node_ids[1:self._tophits_threshold+1]
        self._nodes[nd_id].tophit_ids = set(tophits)

    def _recompute_tophits(self):
        """Recomputes the top-hit candidate set for every active node.
        """
        # for nd_id in self._active_ids:
        #    self._compute_single_tophits_list(nd_id)
    
    def step(self):
        """Executes a single step of the tree-building process.

        In each step the algorithm:
          - Evaluates candidate joining pairs based on the current top-hit sets.
          - Determines the best pair to join (the pair with the smallest distance).
          - Merges the chosen pair into a new node.
          - Periodically refreshes the top-hit candidate sets based on the refresh interval.

        Side Effects:
            Updates internal state including _active_ids, _nodes, _steps, _distance_cache, and _union_find.
        """

        self._steps += 1

        candidate_join_ids = []
        for nd_id1 in self._active_ids:
            best_nd_id = None
            best_distance = float("inf")
            if len(self._nodes[nd_id1].tophit_ids) == 0:
                self._compute_single_tophits_list(nd_id1)
            for nd_id2_temp in self._nodes[nd_id1].tophit_ids:
                nd_id2 = self._union_find.find(nd_id2_temp)
                distance = self._distance_util(nd_id1, nd_id2)
                if distance < best_distance:
                    best_nd_id = nd_id2
                    best_distance = distance
            candidate_join_ids.append((nd_id1, best_nd_id))
        
        # potential postprocessing step with candidate_join_ids

        best_join_ids = min(candidate_join_ids,
                            key=lambda nd_ids: self._distance_util(nd_ids[0], nd_ids[1]))
        nd_id0, nd_id1 = best_join_ids
        self._node_join(nd_id0, nd_id1)

        if self._steps % self._refresh_interval == 0:
            self._recompute_tophits()

    def export_tree(self) -> newick.Node:
        """Exports the constructed tree as a newick.Node with corrected branch distances.

        The tree is represented as an object of type newick.Node. The leaves are labeled with
        the labels of the original sequences, and the edges are weighted according to the distances.

        Returns:
            newick.Node: A newick.Node representing the final phylogenetic tree.
        """

        logging.info("Exporting constructed tree")
        assert len(self._active_ids) == 1
        last_remaining = next(iter(self._active_ids))

        newick_nodes = [newick.Node(i.node_info.label) for i in self._nodes]
        
        def dfs_help(nd_id: NodeID):
            nd = self._nodes[nd_id]
            for child_id in [nd.leftchild_id, nd.rightchild_id]:
                if child_id is None:
                    continue
                raw_dist = self._distance_util(nd_id, child_id)
                distance = CORRECTION(raw_dist)
                dfs_help(child_id)
                newick_nodes[child_id].length = distance
                newick_nodes[nd_id].add_descendant(newick_nodes[child_id])
        dfs_help(last_remaining)
        return newick_nodes[last_remaining]


    def build(self):
        """Executes the full tree-building process.

        Returns:
            List[List[Tuple[int, float]]]: The final adjacency list representation of the
                phylogenetic tree.
        """

        for i in range(self._num_sequences - 1):
            logging.info(f"Step {i+1} of {self._num_sequences-1}")
            self.step()
        return self.export_tree()
