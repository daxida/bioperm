"""Explore permutations.

We shall use the term singlet rather than mononucleotide, doublet in place of
dinucleotide, and triplet instead of trinucleotide. Rather than speak of codons, we
define a triplon as being a member of a set of consecutive nonoverlapping triplets.
These terms are illustrated in figure 1.

References:

Significance of nucleotide sequence alignments...
- https://academic.oup.com/mbe/article/2/6/526/981780
Shuffling biological sequences
- https://www.sciencedirect.com/science/article/pii/S0166218X97814564
Shuffling biological sequences with motif constraints
- https://www.sciencedirect.com/science/article/pii/S1570866707000548
"""

import random
from pprint import pprint
from collections import Counter, defaultdict

S1 = "AGACATAAAGTTCCGTACTGCCGGGAT"
S2 = "AAGTTACGAATACATCCCTGGAGGCGT"  # Doublet-preserving
S3 = "AGTACTGCCGTTCCGGGATAAAGACAT"  # Doublet-and-triplet-preserving
S4 = "AAAGATCCGGTTAGACGGTACTGCCAT"  # Doublet-and-triplon-preserving


E = {
    1: "AA",
    2: "AG",
    11: "AT",
    13: "AC",
    20: "AC",
    23: "AT",
    # "note": "[AA AG]",
    4: "CT",
    9: "CG",
    14: "CC",
    15: "CG",
    21: "CC",
    22: "CA",
    # G
    3: "GC",
    6: "GG",
    7: "GT",
    10: "GA",
    16: "GG",
    17: "GT",
    # "note": "GA",
    # T
    5: "TG",
    8: "TC",
    12: "TA",
    18: "TT",
    19: "TA",
}

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

def is_klet_preserved(s1, s2, k):
    return Counter(chunkify_cons(s1, k)) == Counter(chunkify_cons(s2, k))

def is_klon_preserved(s1, s2, k):
    return Counter(chunkify(s1, k)) == Counter(chunkify(s2, k))

TEST = False

if TEST:
    _SE = [v for _, v in sorted(E.items())]
    SE = _SE[0][0] + "".join(v[1] for v in _SE)
    print(SE, len(SE), len(S1) - 3)
    print(Counter(chunkify_cons(S1, 2)) - Counter(chunkify_cons(SE, 2)))

    print("(S2 vs S1):", is_klet_preserved(S1, S2, 2))
    print("(S3 vs S1):", is_klet_preserved(S1, S3, 2) and is_klet_preserved(S1, S3, 3))
    print("(S4 vs S1):", is_klet_preserved(S1, S4, 2) and is_klon_preserved(S1, S4, 3))


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

def check_if_connected(z_graph: list[str], vertices: list[str]) -> bool:
    """Check if all last edges are connected to sf."""
    assert len(z_graph) == len(vertices) - 1

    uf = UnionFind(len(vertices))
    # print("---")
    # print(z_graph)
    # print(vertices)
    for edge in z_graph:
        fr = edge[:-1]
        to = edge[1:]
        fr_idx = vertices.index(fr)
        to_idx = vertices.index(to)
        # print(fr, to, fr_idx, to_idx)
        uf.union(fr_idx, to_idx)

    # print(uf.group_count())

    return uf.group_count() == 1

def edge_ordering(seq: str, k: int) -> dict[str, list[str]]:
    eord = defaultdict(list)
    for chunk in chunkify_cons(seq, k):
        eord[chunk[:k - 1]].append(chunk)

    # The last vertex must be added if it previously was not.
    eord[seq[-(k-1):]].extend([])

    # The sorting is not needed: for visualization
    eord = dict(sorted(eord.items()))
    return eord


assert len(edge_ordering(S1, 3)) == 16
SEP = "==" * 20


def doublet_preserving_permutation(seq: str, k: int = 2, *, debug: bool = False) -> str:
    # NOTE: k may interfere with keys in dict
    assert k < len(seq)

    # 1. Construct the doublet graph G and edge ordering E(S)
    eord = edge_ordering(seq, k)
    s1 = seq[:k-1]
    sf = seq[-(k-1):]

    if debug:
        print("1.", SEP)
        print(seq)
        pprint(eord)
        print(f"{s1=}, {sf=}")

    # 2. 3. 4.
    last_edges = []
    vertices = list(eord.keys())
    seen = set()

    MAX_ITERATIONS = 100
    while MAX_ITERATIONS > 0:
        last_edges.clear()
        for vertex, edges in eord.items():
            if vertex != sf:
                last_edges.append(random.choice(edges))

        key = tuple(last_edges)
        if key in seen:
            continue
        seen.add(key)

        if check_if_connected(last_edges, vertices):
            break

        MAX_ITERATIONS -= 1

    if not MAX_ITERATIONS:
        original_last_edges = [edges[-1] for edges in eord.values() if edges]
        # We got unlucky and did not even find the original endings once!
        if tuple(original_last_edges) not in seen:
            # msg = f"Exhausted iterations.\nReturning original sequence."
            # print(msg)
            assert False
            return seq
        else:
            msg = f"Exhausted iterations even though we found the original sequence.\n{eord=}\n{seen=}\n{original_last_edges=}\n\nReturning original sequence."
            raise RuntimeError(msg)

    if debug:
        print("5.", SEP)
        print(f"{last_edges=}")

    # Put them at the end
    for last_edge in last_edges:
        vertex = last_edge[:k-1]
        last_edge_idx = eord[vertex].index(last_edge)
        eord[vertex].pop(last_edge_idx)
        eord[vertex].append(last_edge)

    # 5.
    def permute_except_last(lst: list[str]) -> list[str]:
        to_permute = lst[:-1]
        random.shuffle(to_permute)
        return to_permute + [lst[-1]]

    eord_perm = dict()
    for vertex, edges in eord.items():
        if vertex == sf:
            eord_perm[vertex] = edges
            random.shuffle(eord_perm[vertex])
        else:
            perm_edges = permute_except_last(edges)
            eord_perm[vertex] = perm_edges

    if debug:
        print(SEP)
        pprint(eord_perm)

    # 6.
    new_seq = s1
    cur_vertex = s1
    if debug:
        print(k)
        print(cur_vertex, None, new_seq)
    while True:
        try:
            cur_edge = eord_perm[cur_vertex].pop(0)
        except:
            assert cur_vertex == sf
            # print(f"Exiting permutation gracefully: {cur_vertex} was not in edge ordering")
            break
        cur_vertex = cur_edge[1:]
        new_seq += cur_edge[-1]
        if debug:
            print(cur_edge, cur_vertex, new_seq)

    for k, v in eord_perm.items():
        if v:
            msg = f"{k} edge list was not exhausted\n{seq=}\n{last_edges=}"
            raise RuntimeError(msg)
    if len(seq) != len(new_seq):
        msg = f"Expected new_seq to be of len {len(seq)} but was {len(new_seq)}"
        raise RuntimeError(msg)
    if not is_klet_preserved(seq, new_seq, 2):
        msg = "Doublets are not preserved"
        raise RuntimeError(msg)

    return new_seq

# DTP
def doublet_and_triplet_preserving_permutation(seq: str) -> str:
    new_seq = doublet_preserving_permutation(seq, k=3)

    if not is_klet_preserved(seq, new_seq, 2):
        msg = "Doublets are not preserved"
        raise RuntimeError(msg)
    if not is_klet_preserved(seq, new_seq, 3):
        msg = "Triplets are not preserved"
        raise RuntimeError(msg)

    return new_seq

# NUB
def victor(seq: str) -> str:
    new_seq = doublet_preserving_permutation(seq, k=4)

    if not is_klet_preserved(seq, new_seq, 2):
        msg = "Doublets are not preserved"
        raise RuntimeError(msg)
    if not is_klet_preserved(seq, new_seq, 3):
        msg = "Triplets are not preserved"
        raise RuntimeError(msg)
    if not is_klet_preserved(seq, new_seq, 4):
        msg = "Quadruplets are not preserved"
        raise RuntimeError(msg)

    return new_seq

# DtP
def doublet_and_triplon_preserving_permutation(seq: str, *, debug: bool=True) -> str:
    if len(seq) % 3:
        print("Could not DtP permutate, the sequence len is not a multiple of 3")
        return seq

    # 1.
    seq = seq.upper()

    # 2.
    s_star = "".join(seq[i:i+2] + seq[i+2].lower() for i in range(0, len(seq), 3))
    if debug:
        print(f"s ={seq}")
        print(f"s*={s_star}")

    # 3.
    eord = defaultdict(list)
    for chunk in chunkify_cons(s_star, 3):
        eord[chunk[0] + chunk[2]].append(chunk)
    # The sorting is not needed: for visualization
    eord = dict(sorted(eord.items()))
    if debug:
        pprint(eord)

    # 4. Create r
    r = "".join(s_star[i] + s_star[i+2] for i in range(0, len(s_star), 3))

    # 5. permute r
    r_perm = doublet_preserving_permutation(r)
    if debug:
        print(f"r ={r}")
        print(f"r'={r_perm}")

    # 6.
    for _, edges in eord.items():
        random.shuffle(edges)

    new_seq = ""
    for idx in range(0, len(r_perm), 2):
        cur_doublet = r_perm[idx:idx+2]
        try:
            cur_triplon = eord[cur_doublet].pop(0)
        except:
            print(f"Could not find {cur_doublet} in eord")
            break

        new_seq += cur_triplon
        # if debug:
        #     print("> ", r_perm, cur_doublet, "=>", cur_triplon)

    new_seq = new_seq.upper()

    # NOTE: We don't exhaust the edge list as they claim
    # print(eord)
    for k, v in eord.items():
        if v:
            pprint(eord)
            msg = f"{k} edge list was not exhausted"
            raise RuntimeError(msg)

    if len(seq) != len(new_seq):
        msg = f"Expected new_seq to be of len {len(seq)} but was {len(new_seq)}"
        raise RuntimeError(msg)
    if not is_klet_preserved(seq, new_seq, 2):
        msg = "Doublets are not preserved"
        raise RuntimeError(msg)
    if not is_klon_preserved(seq, new_seq, 3):
        msg = "Triplon are not preserved"
        raise RuntimeError(msg)
    # if not is_klet_preserved(seq, new_seq, 3):
    #     print("Ok, DtP did not preserve triplets")

    return new_seq

import string

# For testing
def rand_seq(length: int) -> str:
    return ''.join(random.choices(string.ascii_letters, k=length))

if __name__ == "__main__":
    from Bio import SeqIO

    path = "tests/test_data/test_rand_10000_genomic.fna"
    for record in SeqIO.parse(path, "fasta"):
        # S1 = record.seq

        # print(S1)
        perm = doublet_preserving_permutation(S1)
        # print(perm, "DP")
        perm = doublet_and_triplet_preserving_permutation(S1)
        # print(perm, "DTP")
        perm = doublet_and_triplon_preserving_permutation(S1)
        print(perm, "DtP")
        perm = victor(S1)
        # print(perm, "NUB")

    print("Ok")

