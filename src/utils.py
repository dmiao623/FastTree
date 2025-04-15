import numpy as np

def normalize(A: np.typing.NDArray) -> np.typing.NDArray:
    s = np.sum(A)
    if s == 0:
        return A
    return A / s

class UnionFind:
    def __init__(self, n):
        self._parent = list(range(n))

    def find(self, x):
        if self._parent[x] != x:
            self._parent[x] = self.find(self._parent[x])
        return self._parent[x]

    def union(self, x, y):
        x, y = self.find(x), self.find(y)
        if x == y:
            return
        self._parent[y] = x
