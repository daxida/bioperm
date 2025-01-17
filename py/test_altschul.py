from collections import Counter

from altschul import doublet_preserving_permutation, edge_ordering
from constants import S1, S2, S3, S4
from utils import same_klets, same_klons


def test_altschul_sequences():
    assert same_klets(S1, S2, 2)
    assert same_klets(S1, S3, 2) and same_klets(S1, S3, 3)
    assert same_klets(S1, S4, 2) and same_klons(S1, S4, 3)


def test_edge_ordering():
    assert len(edge_ordering(S1, 3)) == 16


def brute_force_possible_permutations(seq: str, k: int):
    """Get naively all the permutations with same k-lets."""
    import itertools

    all_doublet_perserving_permutations = set()
    for perm in itertools.permutations(seq):
        perm = "".join(perm)
        if same_klets(perm, seq, k):
            all_doublet_perserving_permutations.add(perm)
    print(all_doublet_perserving_permutations)


def test_altschul_klet_uniform():
    """Test that all of possible k-let preserving permutations are found."""
    seq = "ACTAGTAT"

    # Set of all possible permutations found with:
    # brute_force_possible_permutations(seq)
    all_perms = {
        "AGTACTAT",
        "ACTAGTAT",
        "ATAGTACT",
        "ATACTAGT",
        "AGTATACT",
        "ACTATAGT",
    }
    _test_altschul_klet_uniform(seq, all_perms, 2)

    seq = "AAATAAA"
    all_perms = {
        "AAAATAA",
        "AATAAAA",
        "AAATAAA",
    }
    _test_altschul_klet_uniform(seq, all_perms, 3)

    seq = "AAAATAAAA"
    all_perms = {
        "AAAAATAAA",
        "AAATAAAAA",
        "AAAATAAAA",
    }
    _test_altschul_klet_uniform(seq, all_perms, 4)


def _test_altschul_klet_uniform(seq: str, all_perms: set[str], k: int):
    cnt = Counter()
    max_iterations = 5000
    for _ in range(max_iterations):
        res = doublet_preserving_permutation(seq, k)
        cnt[res] += 1
        assert same_klets(seq, res, 1)
        assert same_klets(seq, res, k)

    for perm, amount in cnt.items():
        assert perm in all_perms
        assert abs(amount / max_iterations - 1 / len(all_perms)) < 0.05

    remaining = set(cnt.keys()) - all_perms
    assert len(remaining) == 0
