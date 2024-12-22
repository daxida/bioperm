//! https://www.sciencedirect.com/science/article/pii/S0166218X97814564

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
}
