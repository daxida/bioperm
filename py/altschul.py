"""Explore permutations.

References:

Significance of nucleotide sequence alignments...
- https://academic.oup.com/mbe/article/2/6/526/981780
Shuffling biological sequences
- https://www.sciencedirect.com/science/article/pii/S0166218X97814564
Shuffling biological sequences with motif constraints
- https://www.sciencedirect.com/science/article/pii/S1570866707000548
"""

import random
from collections import defaultdict
from pprint import pprint

from constants import S1
from utils import UnionFind, chunkify_cons, same_klets, same_klons

SEP = "==" * 20


def check_if_connected(z_graph: list[str], vertices: list[str]) -> bool:
    """Check if all last edges are connected to sf."""
    assert len(z_graph) == len(vertices) - 1

    uf = UnionFind(len(vertices))
    for edge in z_graph:
        fr = edge[:-1]
        to = edge[1:]
        fr_idx = vertices.index(fr)
        to_idx = vertices.index(to)
        uf.union(fr_idx, to_idx)

    return uf.group_count() == 1


def edge_ordering(seq: str, k: int) -> dict[str, list[str]]:
    eord = defaultdict(list)
    for chunk in chunkify_cons(seq, k):
        eord[chunk[: k - 1]].append(chunk)

    # The last vertex must be added if it previously was not.
    eord[seq[-(k - 1) :]].extend([])

    # The sorting is not needed: for visualization
    eord = dict(sorted(eord.items()))
    return eord


def doublet_preserving_permutation(seq: str, k: int = 2, *, debug: bool = False) -> str:
    assert k < len(seq)

    # 1. Construct the klet graph G and edge ordering E(S)
    eord = edge_ordering(seq, k)
    s1 = seq[: k - 1]
    sf = seq[-(k - 1) :]

    if debug:
        print("1.", SEP)
        print(seq)
        pprint(eord)
        print(f"{s1=}, {sf=}")

    # 2. For each vertex in G \ s_f, randomly select last edges
    # 3. From the set of last edges, construct z_graph
    #    and assess if connected to s_f
    # 4. If it is not connected, goto (2.)
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

    # Put last edges at the end
    for last_edge in last_edges:
        vertex = last_edge[: k - 1]
        last_edge_idx = eord[vertex].index(last_edge)
        eord[vertex].pop(last_edge_idx)
        eord[vertex].append(last_edge)

    # 5. Randomly permute:
    #    - every edge for s_f
    #    - every edge but the last edge chose for G \ s_f
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

    # 6. Construct the sequence from the permuted edge ordering
    new_seq = s1
    cur_vertex = s1
    if debug:
        print(k)
        print(cur_vertex, None, new_seq)
    while True:
        try:
            cur_edge = eord_perm[cur_vertex].pop(0)
        except IndexError:
            assert cur_vertex == sf
            # print(f"Exiting gracefully: {cur_vertex} was not in edge ordering")
            break
        cur_vertex = cur_edge[1:]
        new_seq += cur_edge[-1]
        if debug:
            print(cur_edge, cur_vertex, new_seq)

    for key, value in eord_perm.items():
        if value:
            msg = f"{key} edge list was not exhausted\n{seq=}\n{last_edges=}"
            raise RuntimeError(msg)
    if len(seq) != len(new_seq):
        msg = f"Expected new_seq to be of len {len(seq)} but was {len(new_seq)}"
        raise RuntimeError(msg)
    if not same_klets(seq, new_seq, 2):
        msg = "Doublets are not preserved"
        raise RuntimeError(msg)

    return new_seq


# DTP
def doublet_and_triplet_preserving_permutation(seq: str) -> str:
    new_seq = doublet_preserving_permutation(seq, k=3)

    if not same_klets(seq, new_seq, 2):
        msg = "Doublets are not preserved"
        raise RuntimeError(msg)
    if not same_klets(seq, new_seq, 3):
        msg = "Triplets are not preserved"
        raise RuntimeError(msg)

    return new_seq


# DtP
def doublet_and_triplon_preserving_permutation(seq: str, *, debug: bool = True) -> str:
    if len(seq) % 3:
        print("Could not DtP permutate, the sequence len is not a multiple of 3")
        return seq

    # 1.
    seq = seq.upper()

    # 2.
    s_star = "".join(seq[i : i + 2] + seq[i + 2].lower() for i in range(0, len(seq), 3))
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
    r = "".join(s_star[i] + s_star[i + 2] for i in range(0, len(s_star), 3))

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
        cur_doublet = r_perm[idx : idx + 2]
        try:
            cur_triplon = eord[cur_doublet].pop(0)
        except IndexError:
            print(f"Could not find {cur_doublet} in eord")
            break

        new_seq += cur_triplon
        # if debug:
        #     print("> ", r_perm, cur_doublet, "=>", cur_triplon)

    new_seq = new_seq.upper()

    # NOTE: We don't exhaust the edge list as they claim
    # print(eord)
    # for k, v in eord.items():
    #     if v:
    #         pprint(eord)
    #         msg = f"{k} edge list was not exhausted"
    #         raise RuntimeError(msg)

    if len(seq) != len(new_seq):
        msg = f"Expected new_seq to be of len {len(seq)} but was {len(new_seq)}"
        raise RuntimeError(msg)
    if not same_klets(seq, new_seq, 2):
        msg = "Doublets are not preserved"
        raise RuntimeError(msg)
    if not same_klons(seq, new_seq, 3):
        msg = "Triplon are not preserved"
        raise RuntimeError(msg)
    # if not is_klet_preserved(seq, new_seq, 3):
    #     print("Ok, DtP did not preserve triplets")

    return new_seq


if __name__ == "__main__":
    # from Bio import SeqIO

    # path = "tests/test_data/test_rand_10000_genomic.fna"
    # for record in SeqIO.parse(path, "fasta"):
    # S1 = record.seq

    # print(S1)
    perm = doublet_preserving_permutation(S1)
    # print(perm, "DP")
    perm = doublet_and_triplet_preserving_permutation(S1)
    # print(perm, "DTP")
    perm = doublet_and_triplon_preserving_permutation(S1)
    print(perm, "DtP")

    print("Ok")
