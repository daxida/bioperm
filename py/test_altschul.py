from constants import S1, S2, S3, S4
from altschul import edge_ordering, doublet_preserving_permutation
from utils import same_klets, same_klons


def test_altschul_sequences():
    assert same_klets(S1, S2, 2)
    assert same_klets(S1, S3, 2) and same_klets(S1, S3, 3)
    assert same_klets(S1, S4, 2) and same_klons(S1, S4, 3)


def test_edge_ordering():
    assert len(edge_ordering(S1, 3)) == 16


def brute_force_possible_permutations(seq: str):
    """Get naively all the permutations with same 2-lets."""
    import itertools

    all_doublet_perserving_permutations = set()
    for perm in itertools.permutations(seq):
        perm = "".join(perm)
        if same_klets(perm, seq, 2):
            all_doublet_perserving_permutations.add(perm)
    print(list(all_doublet_perserving_permutations))


def test_altschul_doublet_preserving_permutation():
    """Test that all of possible 2-let preserving permutations are found."""
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

    for _ in range(1000):
        res = doublet_preserving_permutation(seq, k=2)
        all_perms.discard(res)
        assert same_klets(seq, res, 2)

    assert len(all_perms) == 0
