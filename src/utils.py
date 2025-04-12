import numpy as np

def normalize(A: np.typing.NDArray) -> np.typing.NDArray:
    s = np.sum(A)
    if s == 0:
        return A
    return A / s
