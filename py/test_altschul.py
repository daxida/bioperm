from constants import S1, S2, S3, S4
from altschul import edge_ordering
from utils import same_klets, same_klons


def test_altschul_sequences():
    assert same_klets(S1, S2, 2)
    assert same_klets(S1, S3, 2) and same_klets(S1, S3, 3)
    assert same_klets(S1, S4, 2) and same_klons(S1, S4, 3)


def test_edge_ordering():
    assert len(edge_ordering(S1, 3)) == 16
