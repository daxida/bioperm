"""Swap and euler algorithms.

Shuffling biological sequences
- https://www.sciencedirect.com/science/article/pii/S0166218X97814564
"""

# from constants import S1
import random
from utils import same_klets
from altschul import edge_ordering
from pprint import pprint


def is_k_cyclic(seq: str, k: int) -> bool:
    """Does no check if k makes sense."""
    return seq[: (k - 1)] == seq[-(k - 1) :]


def random_rotation(seq: str, k: int, *, m: int | None = None) -> str:
    """Preserves k-lets.

    Note that this is not a proper rotation. It does not
    conserve the characters of the sequence.

    m can be passed as an int for testing.
    """
    assert is_k_cyclic(seq, k)

    n = len(seq)
    if m is None:
        m = random.randint(k, n)
    assert k <= m <= n

    rotated = seq[m - 1 :] + seq[k - 1 : m]
    for i in range(m + 1, m + k - 1):
        idx = i if i <= n else i % n + k - 1
        rotated += seq[idx - 1]

    assert len(seq) == len(rotated)
    assert same_klets(seq, rotated, k)

    return rotated


def markov_transition(seq: str, k: int) -> str | None:
    """Return the swapped sequence if the checks are verified, or else None."""
    n = len(seq)
    positions = random.sample(range(0, n - k + 3), 4)
    assert all(0 <= x <= n - k + 2 for x in positions)
    positions.sort()  # Ensure a < b < c < d
    a, b, c, d = positions

    ss1 = seq[a : a + k - 1]
    ss2 = seq[b : b + k - 1]
    ss3 = seq[c : c + k - 1]
    ss4 = seq[d : d + k - 1]
    assert len(ss1) == k - 1

    if ss1 == ss3 and ss2 == ss4:
        # print(a, b, c, d)
        overlap = b + k - 1 >= c
        res = (
            seq[:a]
            + seq[c : d + k - 1]  # swapped
            + seq[b + k - 1 : c]
            # + seq[a : min(c, b + k - 1)]  # swapped
            + seq[a : b + k - 1]  # swapped
            + seq[d + k - 1 :]
        )
        # TODO: finish the overlap case (seq AAATAAA).
        diagnostic = (
            f"{seq=} {len(seq)=} != {len(res)=} {res=}\n"
            f"{overlap=} {b + k - 1} >= {c}\n"
            f"{a} {b} {c} {d}\n"
            f"[EXPECTED] Swapping {seq[a : b + k - 1]} with {seq[c : d + k - 1]}\n"
            # Our implementation
            f"{ss1} [{a}..{a + k - 1}] == {ss3} [{c}..{c + k - 1}]\n"
            f"{ss2} [{b}..{b + k - 1}] == {ss4} [{d}..{d + k - 1}]\n"
            # They mention this invariant in the article
            f"{seq[c : b + k - 1]} [{c}..{b + k - 1}]"
        )
        assert len(seq) == len(res), diagnostic
        return res
    else:
        return None


def swap_algorithm(seq: str, k: int) -> str:
    """Section 3."""
    assert k < len(seq)
    was_cyclic = is_k_cyclic(seq, k)
    if was_cyclic:
        seq = random_rotation(seq, k)
        # raise NotImplementedError(f"{k}-cyclic {seq=}")

    new_seq = None
    max_iterations = 5000
    for _ in range(max_iterations):
        new_seq = markov_transition(seq, k)
        if new_seq is not None:
            break
    # Max iterations reached
    if new_seq is None:
        # Cf. AATAA -> rotated to -> TAAAT (with no possible transition)
        return seq
        # raise RuntimeError(f"{max_iterations=} reached in swap_algorithm for {seq=}")

    if was_cyclic:
        assert is_k_cyclic(seq, k)
    assert len(new_seq) == len(seq), (
        f"{seq=} length != {new_seq=} length: {len(seq)} != {len(new_seq)}"
    )
    assert same_klets(seq, new_seq, k)

    return new_seq


def euler_algorithm(seq: str, k: int) -> str:
    """Euler algorithm for generating a uniform random permutation of seq.

    1. Construct the digraph D_k(seq) from klet counts.
    2. If seq is cyclic...
    """
    assert k < len(seq)

    # 1. Construct D_k
    # eord = defaultdict(list)
    # for chunk in chunkify_cons(seq, k):
    #     eord[chunk[: k - 1]].append(chunk[1:])
    eord = edge_ordering(seq, k)

    # 2. If acyclic, add dummy edge.
    #    If not acyclic, random rotation.
    if is_k_cyclic(seq, k):
        raise NotImplementedError

    fst = seq[: k - 1]
    lst = seq[-(k - 1) :]
    dummy_vertex = f"X{fst}"
    eord[lst].append(dummy_vertex)
    pprint(eord)

    # 3. Random walk. Add first-encountered edges to T
    cur_vertex = lst
    t = [cur_vertex]
    seen = {cur_vertex}
    while len(seen) < len(eord):
        cur_edge = eord[cur_vertex][0]
        cur_vertex = cur_edge[1:]
        seen.add(cur_vertex)
        t.append(cur_vertex)

    print(seq, t)
    # 4. Randomly order all arcs not in T

    # 5. Read desired sequence starting at fst

    # Assumes we reordered the edges..
    eord_perm = eord
    eord_perm[lst].remove(dummy_vertex)
    new_seq = fst
    cur_vertex = fst
    while True:
        try:
            cur_edge = eord_perm[cur_vertex].pop(0)
        except IndexError:
            msg = f"Last vertex {cur_vertex} should be equal to {lst}"
            assert cur_vertex == lst, msg
            # print(f"Exiting gracefully: {cur_vertex} was not in edge ordering")
            break
        cur_vertex = cur_edge[1:]
        new_seq += cur_edge[-1]
        # print(cur_edge, cur_vertex, new_seq)

    print(f"{new_seq=}")
    assert same_klets(seq, new_seq, k)

    return new_seq


def main():
    S1 = "ATCAGCAAC"
    S1 = "AAATAAA"
    print(S1)
    new_seq = swap_algorithm(S1, 2)
    print(new_seq)


if __name__ == "__main__":
    main()
