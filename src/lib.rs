/// Chunkify function for non-overlapping blocks
pub fn chunkify(seq: &str, chunk_size: usize) -> Vec<String> {
    (0..seq.len())
        .step_by(chunk_size)
        .map(|i| seq[i..std::cmp::min(i + chunk_size, seq.len())].to_string())
        .collect()
}

// Chunkify function for overlapping blocks
pub fn chunkify_cons(seq: &str, chunk_size: usize) -> Vec<String> {
    if seq.len() < chunk_size {
        return Vec::new();
    }
    (0..=seq.len() - chunk_size)
        .map(|i| seq[i..i + chunk_size].to_string())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chunkify() {
        assert_eq!(chunkify("ABCDEFG", 3), vec!["ABC", "DEF", "G"]);
        assert_eq!(chunkify("ABCDE", 2), vec!["AB", "CD", "E"]);
        assert_eq!(chunkify("ABCDE", 5), vec!["ABCDE"]);
        assert_eq!(chunkify("A", 1), vec!["A"]);
        assert_eq!(chunkify("ABCDE", 1), vec!["A", "B", "C", "D", "E"]);
    }

    #[test]
    fn test_chunkify_cons() {
        assert_eq!(chunkify_cons("ABCDE", 2), vec!["AB", "BC", "CD", "DE"]);
        assert_eq!(chunkify_cons("ABCDE", 3), vec!["ABC", "BCD", "CDE"]);
        assert_eq!(chunkify_cons("A", 1), vec!["A"]);
    }
}
