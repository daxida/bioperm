//! https://www.sciencedirect.com/science/article/pii/S0166218X97814564
use crate::utils::same_klets;
use rand::Rng;

fn is_k_cyclic(seq: &str, k: usize) -> bool {
    let prefix = &seq[..k - 1];
    let suffix = &seq[seq.len() - (k - 1)..];
    prefix == suffix
}

fn random_rotation(seq: &str, k: usize, m: Option<usize>) -> String {
    assert!(is_k_cyclic(seq, k), "Sequence is not k-cyclic");
    let n = seq.len();

    let m = m.unwrap_or_else(|| rand::thread_rng().gen_range(k..=n));
    assert!(k <= m && m <= n);

    let mut rotated = String::new();
    rotated.push_str(&seq[m - 1..]);
    rotated.push_str(&seq[k - 1..m]);

    for i in (m + 1)..(m + k - 1) {
        let idx = if i <= n { i } else { i % n + k - 1 };
        rotated.push(seq.chars().nth(idx - 1).unwrap());
    }

    assert_eq!(seq.len(), rotated.len());
    assert!(same_klets(seq, &rotated, k));

    rotated
}

#[allow(dead_code)]
fn kandel(seq: &str) -> &str {
    // Implement the exact algorithm p.181
    seq
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::same_klets;

    #[test]
    fn test_same_klets_two() {
        let chunk_size = 2;
        let s1 = "AGACATAAAGTTCCGTACTGCCGGGAT";
        let s2 = kandel(s1);
        assert!(same_klets(s1, &s2, chunk_size));
    }

    #[test]
    fn test_is_k_cyclic() {
        assert!(is_k_cyclic("AATAA", 2));
        assert!(is_k_cyclic("AATAA", 3));
    }

    fn test_random_rotation(seq: &str, k: usize, expected: &[&str]) {
        let n = seq.len();
        let mut expected_seqs: Vec<_> = expected.iter().map(|s| s.to_string()).collect();
        let mut received_seqs: Vec<_> = (k..=n).map(|m| random_rotation(seq, k, Some(m))).collect();

        expected_seqs.sort();
        received_seqs.sort();

        assert_eq!(expected_seqs.len(), n - k + 1);
        assert_eq!(expected_seqs, received_seqs);
    }

    #[test]
    fn test_random_rotation_base() {
        let seq = "AATAA";
        let k = 2;
        let expected = ["ATAAA", "TAAAT", "AAATA", "AATAA"];
        test_random_rotation(seq, k, &expected);
    }

    #[test]
    fn test_random_rotation_base_two() {
        let seq = "ACGTAC";
        let k = 3;
        let expected = ["ACGTAC", "CGTACG", "GTACGT", "TACGTA"];
        test_random_rotation(seq, k, &expected);
    }
}
