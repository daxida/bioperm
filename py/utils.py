import random
import string
from collections import Counter


class UnionFind:
    def __init__(self, n):
        self.n = n
        self.parents = [-1] * n
        self.group = n

    def find(self, x):
        if self.parents[x] < 0:
            return x
        else:
            self.parents[x] = self.find(self.parents[x])
            return self.parents[x]

    def union(self, x, y):
        x = self.find(x)
        y = self.find(y)

        if x == y:
            return
        self.group -= 1
        if self.parents[x] > self.parents[y]:
            x, y = y, x

        self.parents[x] += self.parents[y]
        self.parents[y] = x

    def group_count(self):
        return self.group


def rand_seq(length: int) -> str:
    """Testing function."""
    return "".join(random.choices(string.ascii_letters, k=length))


def chunkify(seq: str, chunk_size: int) -> list[str]:
    """Chunkify consecutive, non overlapping blocks.

    chunkify("ABCDEFG", 3) → ABC DEF G

    Identical to https://docs.python.org/3/library/itertools.html#itertools.batched
    """
    return [seq[fr : fr + chunk_size] for fr in range(0, len(seq), chunk_size)]


def chunkify_cons(seq: str, chunk_size: int) -> list[str]:
    """Chunkify consecutive, overlapping blocks.

    chunkify_cons("ABCDE", 2) → AB BC CD DE
    chunkify_cons("ABCDE", 3) → ABC BCD CDE

    Similar to ruby "each_cons".
    If chunk_size = 2, then it is identical to
    https://docs.python.org/3/library/itertools.html#itertools.pairwise
    """
    if len(seq) < chunk_size:
        return []
    return [seq[fr : fr + chunk_size] for fr in range(len(seq) - chunk_size + 1)]


def same_klets(s1, s2, k):
    """True iif s1 and s2 have same klets."""
    return Counter(chunkify_cons(s1, k)) == Counter(chunkify_cons(s2, k))


def same_klons(s1, s2, k):
    """True iif s1 and s2 have same klons."""
    return Counter(chunkify(s1, k)) == Counter(chunkify(s2, k))
