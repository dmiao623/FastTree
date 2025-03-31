from typing import List, Tuple, TypeVar

from sequence import AbstractSequence
Sequence = TypeVar("Sequence", bound=AbstractSequence)

class Alignment:
    def __init__(self, alignment: List[Sequence]):
        if not alignment:
            raise ValueError("Profile must be initialized with at least one sequence.")

        self._alignment_size = len(alignment)

        # [TODO]: compute multi-alignment

        self._alignment = None
        self._alignment_matrix = self._compute_alignment_matrix()
        self._alignment_length = len(self._alignment_matrix[0])

    def _compute_alignment_matrix(self):
        return [[]]

    @property
    def alignment(self): return self._alignment

    @property
    def alignment_size(self): return self._alignment_size

    @property
    def alignment_length(self): return self._alignment_length

    def __getitem__(self, key: Tuple[int, int]):
        sequence_id = key[0]
        position_id = key[1]
        return self._alignment_matrix[sequence_id][position_id]
