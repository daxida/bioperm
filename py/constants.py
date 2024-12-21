"""Altschul' article sequeces."""

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

# _SE = [v for _, v in sorted(E.items())]
# SE = _SE[0][0] + "".join(v[1] for v in _SE)
# print(SE, len(SE), len(S1) - 3)
# print(Counter(chunkify_cons(S1, 2)) - Counter(chunkify_cons(SE, 2)))
