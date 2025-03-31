import numpy as np

from typing import Union
from numpy.typing import NDArray

from constants import DNAConstants

class Sequence:
    """Class representing a generic biological sequence.

    Attributes:
        _sequence (NDArray[int]): 1D array of encoded sequence values.
        _sequence_length (int): length of sequence, including gaps.

    Properties:
        sequence (NDArray[int]): the stored sequence as a 1D array of integers
        sequence_length (int): the length of the stored sequence
    """


    def __init__(self, sequence: Union[NDArray[int], str]): 
        """ Constructs an AbstractSequence object.

        Args:
            sequence (NDArray[int]): 1D Numpy array of sequence characters, represented as integers
                in [0, k) where k is the size of the alphabet.

        Raises:
                ValueError: Raised if sequence does not have dimension 1 or is empty.
        """

        if isinstance(sequence, str):
            try: 
                sequence = np.array([DNAConstants.NUCLEOTIDES.index(c) for c in sequence])
            except ValueError:
                raise ValueError("Provided sequence contains invalid nucleotide characters.")

        if arr.ndim != 1:
            raise ValueError("Invalid dimension of sequence provided.")
        if arr.size == 1:
            raise ValueError("Provided sequence is empty.")

        self._sequence_length = len(sequence)
        self._sequence = sequence



    @property
    def sequence_length(self) -> int: return self._sequence_length

    @property
    def sequence(self) -> NDArray[int]: return self._sequence

    def __str__(self) -> str:
        return "".join([constants.ALPHABET[c] for c in self._sequence])

    def __iter__(self): return iter(self._sequence)
