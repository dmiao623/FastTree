import numpy as np

def _normalize(A: np.typing.NDArray) -> np.typing.NDArray:
    s = np.sum(A)
    if s == 0:
        return A
    return A / s
