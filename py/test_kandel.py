from kandel import random_rotation_m, is_k_cyclic
from utils import is_klet_preserved


def test_random_rotation() -> None:
    seq = "ACGTAC"
    k = 3
    for m in range(k, len(seq) + 1):
        rotated = random_rotation_m(seq, k, m)
        assert is_k_cyclic(rotated, k)
        assert is_klet_preserved(seq, rotated, k)


def test_random_rotation_different_k() -> None:
    seq = "ACGTACG"
    k = 4
    for m in range(k, len(seq) + 1):
        rotated = random_rotation_m(seq, k, m)
        assert is_k_cyclic(rotated, k)
        assert is_klet_preserved(seq, rotated, k)
