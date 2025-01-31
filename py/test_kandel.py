from collections import Counter

from kandel import is_k_cyclic, random_rotation, swap_algorithm
from utils import same_klets


def test_random_rotation_base() -> None:
    seq = "AATAA"  # 2 and 3-cyclic
    #      12345
    k = 2
    n = len(seq)
    # m is in k..=n (2..=5)
    expected_seqs = {
        "seq_m2": "ATAA" + "A",
        "seq_m3": "TAA" + "AT",
        "seq_m4": "AA" + "ATA",
        "seq_m5": "A" + "ATAA",
    }.values()
    received_seqs = [random_rotation(seq, k, m=m) for m in range(k, n + 1)]
    assert {*expected_seqs} == {*received_seqs}


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


def test_klet_uniform():
    """Test that all of possible k-let preserving permutations are found."""
    seq = "ACTAGTAT"

    all_perms = {
        "AGTACTAT",
        "ACTAGTAT",
        "ATAGTACT",
        "ATACTAGT",
        "AGTATACT",
        "ACTATAGT",
    }
    _test_klet_uniform(seq, all_perms, 2)


def _test_klet_uniform(seq: str, all_perms: set[str], k: int):
    cnt = Counter()
    max_iterations = 1000
    for _ in range(max_iterations):
        res = swap_algorithm(seq, k)
        cnt[res] += 1
        if not is_k_cyclic(seq, k):
            assert same_klets(seq, res, 1)
        assert same_klets(seq, res, k)

    for perm, amount in cnt.items():
        assert perm in all_perms, perm
        assert abs(amount / max_iterations - 1 / len(all_perms)) < 0.05

    remaining = set(cnt.keys()) - all_perms
    assert len(remaining) == 0
