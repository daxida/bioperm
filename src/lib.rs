use std::collections::HashMap;
use std::hash::Hash;

/// Chunkify function for non-overlapping blocks
///
/// "ABCDEFG", 3 > ("ABC", "DEF", "G")
pub fn each_step(seq: &str, step: usize) -> impl Iterator<Item = &str> {
    (0..seq.len())
        .step_by(step)
        .map(move |i| &seq[i..std::cmp::min(i + step, seq.len())])
}

/// Chunkify function for overlapping blocks
///
/// "ABCDEFG", 3 > ("ABC", "BCD", "CDE", "DEF", "EFG")
pub fn each_cons(seq: &str, step: usize) -> impl Iterator<Item = &str> {
    (0..=seq.len() - step).map(move |i| &seq[i..i + step])
}

fn same_frequencies<I, T>(it1: I, it2: I) -> bool
where
    I: Iterator<Item = T>,
    T: Eq + Hash,
{
    it1.zip(it2)
        .fold(HashMap::new(), |mut cnt, (item1, item2)| {
            *cnt.entry(item1).or_insert(0) += 1;
            *cnt.entry(item2).or_insert(0) -= 1;
            cnt
        })
        .values()
        .all(|&v| v == 0)
}

pub fn same_klets(s1: &str, s2: &str, k: usize) -> bool {
    s1.len() == s2.len() && same_frequencies(each_cons(s1, k), each_cons(s2, k))
}

pub fn same_klons(s1: &str, s2: &str, k: usize) -> bool {
    s1.len() == s2.len() && same_frequencies(each_step(s1, k), each_step(s2, k))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chunkify() {
        assert_eq!(
            each_step("ABCDEFG", 3).collect::<Vec<_>>(),
            vec!["ABC", "DEF", "G"]
        );
        assert_eq!(
            each_step("ABCDE", 2).collect::<Vec<_>>(),
            vec!["AB", "CD", "E"]
        );
    }

    #[test]
    fn test_chunkify_cons() {
        assert_eq!(
            each_cons("ABCDE", 2).collect::<Vec<_>>(),
            vec!["AB", "BC", "CD", "DE"]
        );
        assert_eq!(
            each_cons("ABCDE", 3).collect::<Vec<_>>(),
            vec!["ABC", "BCD", "CDE"]
        );
    }

    #[test]
    fn test_same_klets() {
        let s1 = "AGACATAAAGTTCCGTACTGCCGGGAT";
        let s4 = "AAAGATCCGGTTAGACGGTACTGCCAT";
        assert!(same_klets(s1, s4, 2));
    }

    #[test]
    fn test_same_klons() {
        let s1 = "AGACATAAAGTTCCGTACTGCCGGGAT";
        let s4 = "AAAGATCCGGTTAGACGGTACTGCCAT";
        assert!(same_klons(s1, s4, 3));
    }
}
