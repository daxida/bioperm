from constants import S1, S2, S3, S4
from main import edge_ordering
from utils import is_klet_preserved, is_klon_preserved


def test_altschul_sequences():
    assert is_klet_preserved(S1, S2, 2)
    assert is_klet_preserved(S1, S3, 2) and is_klet_preserved(S1, S3, 3)
    assert is_klet_preserved(S1, S4, 2) and is_klon_preserved(S1, S4, 3)


def test_edge_ordering():
    assert len(edge_ordering(S1, 3)) == 16
