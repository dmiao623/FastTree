from typing import Dict

from profile import Profile

class Alignment:
    """
    An alignment of DNA sequences.

    [ToDo]: is this class even necessary?

    Properties:
        alignment (Dict[str, str]): a dictionary with the multiple alignment;
            the labels of the sequences are the keys, the sequences are the values
        alignment_size (int): how many sequences are aligned
        alignment_length (int): the length of the aligned sequences
        profile_dict (Dict[str, Profile]): a dictionary with the profile for each sequence
    """
    def __init__(self, alignment: Dict[str, str]):
        """
        The input is given as a dictionary with the labels as the keys and
        the sequences as the values.
        """
        if not alignment:
            raise ValueError("Alignment must be initialized with at least one sequence.")

        self._alignment = alignment
        self._alignment_size = len(alignment)
        # The line below gets any string in the dict and takes its length.
        self._alignment_length = len(next(iter(alignment.values())))

        if not all(len(seq) == self._alignment_length for seq in alignment.values()):
            raise ValueError("Sequences in alignment do not all have the same length.")

        self._profile_dict = {label: Profile.from_aligned_sequence(seq) for label, seq in alignment.items()}

    @property
    def alignment(self): return self._alignment

    @property
    def alignment_size(self): return self._alignment_size

    @property
    def alignment_length(self): return self._alignment_length

    @property
    def profile_dict(self): return self._profile_dict
