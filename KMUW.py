'''Chatgpt code for the article 
Shuffling biological sequences
from D. Kandel”, Y. MatiasbT*, R. Unger”, P. Winklerb'''
import random

def find_flanked_substrings(sequence, k):
    """
    Identify all pairs of substrings that are disjoint, flanked by the same (k-1)-lets,
    and satisfy the swapping conditions.
    """
    n = len(sequence)
    candidates = []

    for a in range(n - 2 * k + 2):
        for b in range(a + k - 1, n - k + 1):
            left_flank = sequence[a:a + k - 1]
            right_flank = sequence[b:b + k - 1]

            for c in range(b + k - 1, n - k + 1):
                for d in range(c + k - 1, n):
                    if (sequence[a:a + k - 1] == sequence[c:c + k - 1] and
                            sequence[b:b + k - 1] == sequence[d:d + k - 1]):
                        candidates.append(((a, b), (c, d)))

    return candidates

def swap_substrings(sequence, indices):
    """
    Perform the swap of substrings defined by indices.
    """
    (a, b), (c, d) = indices
    substring_1 = sequence[a:b + 1]
    substring_2 = sequence[c:d + 1]

    new_sequence = (
        sequence[:a] + substring_2 + sequence[b + 1:c] + substring_1 + sequence[d + 1:]
    )

    return new_sequence

def random_rotation(sequence, k):
    """
    Perform a random rotation for cyclic sequences.
    """
    n = len(sequence)
    m = random.randint(k, n)
    rotated_sequence = sequence[m:] + sequence[:m]
    return rotated_sequence

def swap_algorithm(sequence, k, max_iterations=1000):
    """
    Run the swapping algorithm on the given sequence to preserve k-let counts.
    """
    n = len(sequence)

    for _ in range(max_iterations):
        # Randomly choose positions a, b, c, d
        a, b, c, d = sorted(random.sample(range(n - k + 2), 4))
        if a < b < c < d:
            if (sequence[a:a + k - 1] == sequence[c:c + k - 1] and
                    sequence[b:b + k - 1] == sequence[d:d + k - 1]):
                sequence = swap_substrings(sequence, ((a, b), (c, d)))

    return sequence

# Example usage
if __name__ == "__main__":
    original_sequence = "TTACACTGATTCAAGTTAAT"
    k = 2

    # Random rotation if the sequence is cyclic
    if original_sequence[:k - 1] == original_sequence[-(k - 1):]:
        original_sequence = random_rotation(original_sequence, k)

    print("Original Sequence:", original_sequence)
    swapped_sequence = swap_algorithm(original_sequence, k)
    print("Swapped Sequence:", swapped_sequence)
