import abc
import numpy as np

from typing import List, Union
from numpy.typing import NDArray

from constants import DNAConstants

class AbstractSequence(abc.ABC):
    """
    Abstract base class representing a generic biological sequence.

    Attributes:
        _sequence (NDArray[int]): 1D array of encoded sequence values.

    Properties:
        sequence (NDArray[int]): The stored sequence as a 1D array of integers.
        alphabet (NDArray[np.str_]): The alphabet corresponding to the sequence (abstract).
    """


    def __init__(self, sequence: NDArray[int]): 
        """ Constructs an AbstractSequence object.

        Args:
            sequence (NDArray[int]): 1D Numpy array of sequence characters, represented as integers
                in [0, k) where k is the size of the alphabet.

        Raises:
                ValueError: Raised if sequence does not have dimension 1 or is empty.
        """

        if arr.ndim != 1:
            raise ValueError("Invalid dimension of sequence provided.")
        if arr.size == 1:
            raise ValueError("Provided sequence is empty.")

        self._sequence = sequence

    @property
    def sequence(self) -> NDArray[int]: return self._sequence

    @property
    @abc.abstractmethod
    def alphabet(self) -> NDArray[np.str_]: pass

    @abc.abstractmethod
    def __str__(self) -> str:
        pass

class DNASequence(AbstractSequence):
    """
    Concrete implementation of AbstractSequence for DNA sequences.

    Attributes:
        _sequence (NDArray[int]): Encoded DNA sequence.

    Properties:
        alphabet (NDArray[np.str_]): Nucleotide alphabet ["A", "C", "G", "T"].
    """

    def __init__(self, sequence: Union[str, NDArray[int]]):
        """ Constructs a DNDASequence object.
        
        Args:
            sequence (NDArray[int] | str): 1D Numpy array of sequence characters, represented a string 
                of "AGCT" or as an array of integers between [0, 3]. 

        Raises:
            ValueError: Raised if provided sequence is invalid.
        """
        
        if isinstance(sequence, str):
            try: 
                sequence = np.array([DNAConstants.NUCLEOTIDES.index(c) for c in sequence])
            except ValueError:
                raise ValueError("Provided sequence contains invalid nucleotide characters.")

        super().__init__(sequence)

    @property
    def alphabet(self) -> NDarray[np.str_]: return DNAConstants.NUCLEOTIDES

    def __str__(self) -> str:
        return "".join([DNAConstants.NUCLEOTIDES[c] for c in self._sequence])
