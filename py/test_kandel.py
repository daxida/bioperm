from kandel import random_rotation, is_k_cyclic
from utils import same_klets


def test_random_rotation() -> None:
    seq = "ACGTAC"
    k = 3
    for m in range(k, len(seq) + 1):
        rotated = random_rotation(seq, k, m=m)
        assert is_k_cyclic(rotated, k)
        assert same_klets(seq, rotated, k)


def test_random_rotation_different_k() -> None:
    seq = "ACGTACG"
    k = 4
    for m in range(k, len(seq) + 1):
        rotated = random_rotation(seq, k, m=m)
        assert is_k_cyclic(rotated, k)
        assert same_klets(seq, rotated, k)
