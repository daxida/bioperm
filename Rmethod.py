from itertools import product
import random
ALPHABET = "ATGC"

'''Old in house method for k-mer conservation'''
def method_r(
    seq: str,
    chunk_size: int,
    *,
    # Whether to keep looping over all split patterns.
    # This being true guarantees more randomness.
    loop: bool = True,
) -> str:
    """Guarantee that all possible consecutive overlapping k-mer are conserved with an specific k-mer.

    > CTATTGGCGTCCACCATTCCTTCGATTATCGCGCCCACTC

    > cut it into n chunks where an specific (k-1)-mer it is found (AT)
    >   CT AT  TGGCGTCCACC  AT TCCTTCG  AT T  AT CGCGCCCACTC
    > ['CT',  'TGGCGTCCACC',  'TCCTTCG',  'T',  'CGCGCCCACTC']

    > find min(x,(n-2)!-1) different configuration without permuting first
    > and final chunk, here (4-2)!-1 = 2!-1 = 1
    > where x it is the number of different configurations we want.
    > ['CT',  'TCCTTCG', 'TGGCGTCCACC', 'T', 'CGCGCCCACTC']
    >   CT  AT TCCTTCG AT TGGCGTCCACC AT T AT CGCGCCCACTC
    """
    # Minimun chunk for the algo to work
    min_chunk = 3

    split_patterns = ["".join(kmer) for kmer in product(ALPHABET, repeat=chunk_size - 1)]
    random.shuffle(split_patterns)
    permutation = "".join(seq)

    for pattern in split_patterns:
        seq = permutation

        chunks_new = seq.split(pattern)

        if len(chunks_new) < min_chunk:
            continue

        fst, chunks_shuffle, lst = chunks_new[0], chunks_new[1:-1], chunks_new[-1]
        random.shuffle(chunks_shuffle)

        permutation = pattern.join(chunks_shuffle)
        # add A at the end to chunks[1:], T at the begining to chunks[:-1]
        permutation = pattern.join([fst, permutation, lst])
        if not loop:
            return permutation

    return permutation


if __name__ == "__main__":
    original_sequence = "TTACACTGATTCAAGTTAAT"
    k = 2 

    print("Original Sequence:", original_sequence)
    swapped_sequence = method_r(original_sequence, k)
    print("Swapped Sequence:", swapped_sequence)