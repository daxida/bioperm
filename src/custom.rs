use itertools::{Itertools, MultiProduct};
use rand::seq::SliceRandom;

fn product_repeat<I>(it: I, repeat: usize) -> MultiProduct<I>
where
    I: Iterator + Clone,
    I::Item: Clone,
{
    std::iter::repeat(it).take(repeat).multi_cartesian_product()
}

/// Custom permutation. Preserves k-mers.
///
/// ## Example
///
/// `CTATTTCACCATTCGATT` with k=3
///
/// Split by a specific (k-1)-mer. Let's say AT:
/// `  CT AT  TTCACC  AT TCG  AT T`
/// `['CT',  'TTCACC',  'TCG',  'T']`
///
/// Then swap some of the chunks, and join by AT:
/// `['CT',  'TCG',  'TTCACC',  'T']`
/// `  CT AT  TCG  AT TTCACC  AT T`
pub fn method_r(seq: &str, k: usize) -> String {
    let mut rng = rand::thread_rng();

    for split in product_repeat("ACGT".chars(), k - 1) {
        let split_str: String = split.into_iter().collect();
        let chunks = seq.split(&split_str).collect::<Vec<_>>();

        if chunks.len() > 2 {
            let fst = chunks[0];
            let lst = chunks[chunks.len() - 1];
            let mut chunks_shuffle = chunks[1..chunks.len() - 1].to_vec();

            chunks_shuffle.shuffle(&mut rng);

            let perm = [fst, &chunks_shuffle.join(&split_str), lst];
            return perm.join(&split_str);
        }
    }

    unreachable!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::same_klets;

    #[test]
    fn test_same_klets_two() {
        let chunk_size = 2;
        let s1 = "AGACATAAAGTTCCGTACTGCCGGGAT";
        let s2 = method_r(s1, chunk_size);
        assert!(same_klets(s1, &s2, chunk_size));
    }

    #[test]
    fn test_same_klets_three() {
        let chunk_size = 3;
        let s1 = "AGACATAAAGTTCCGTACTGCCGGGAT";
        let s2 = method_r(s1, chunk_size);
        assert!(same_klets(s1, &s2, chunk_size));
    }
}
