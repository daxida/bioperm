//! https://www.sciencedirect.com/science/article/pii/S0166218X97814564
use crate::utils::same_klets;
use rand::{seq::SliceRandom, Rng};

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

fn markov_transition(seq: &str, k: usize) -> Option<String> {
    use rand::thread_rng;
    let n = seq.len();
    if n < k {
        return None;
    }
    let mut rng = thread_rng();
    let mut positions: Vec<usize> = (0..n - k + 3).collect();
    positions.shuffle(&mut rng);
    positions.truncate(4);
    positions.sort();
    
    let (a, b, c, d) = (positions[0], positions[1], positions[2], positions[3]);
    
    let ss1 = &seq[a..a + k - 1];
    let ss2 = &seq[b..b + k - 1];
    let ss3 = &seq[c..c + k - 1];
    let ss4 = &seq[d..d + k - 1];
    assert_eq!(ss1.len() == (k-1));
    
    if ss1 == ss3 && ss2 == ss4 {
        let mut res = String::new();
        res.push_str(&seq[..a]);
        res.push_str(&seq[c..d + k - 1]);
        res.push_str(&seq[b + k - 1..c]);
        res.push_str(&seq[a..b + k - 1]);
        res.push_str(&seq[d + k - 1..]);
        
        assert_eq!(seq.len(), res.len());
        return Some(res);
    }
    None
}

fn swap_algorithm(seq: &str, k: usize) -> String {
    assert!(k < seq.len());
    let was_cyclic: bool = is_k_cyclic(seq, k);

    let mut new_seq = None;

    if (was_cyclic){
        new_seq = random_rotation(seq, k, kom)e(m)
    }

    let max_iterations = 5000;
    for _ in 0..max_iterations {
        if let Some(res) = markov_transition(seq, k) {
            new_seq = Some(res);
            break;
        }
    }  
    
    if was_cyclic {
        assert!(is_k_cyclic(&new_seq, k));
    }
    
    let new_seq = new_seq.unwrap_or_else(|| seq.to_string());
    
    assert_eq!(seq.len(), new_seq.len(), "Length mismatch between original and transformed sequence");
    assert!(same_klets(seq, &new_seq, k), "K-let mismatch between original and transformed sequence");
    
    new_seq
}

// fn euler_algorithm(seq: &str, k: usize) -> String {
//     assert!(k < seq.len());
//     let mut eord = edge_ordering(seq, k);
//     let fst = &seq[..k - 1];
//     let lst = &seq[seq.len() - k + 1..];
//     let dummy_vertex = format!("X{}", fst);
//     eord.entry(lst.to_string()).or_default().push(dummy_vertex.clone());
    
//     let mut cur_vertex = lst.to_string();
//     let mut t = vec![cur_vertex.clone()];
//     let mut seen: HashSet<String> = HashSet::new();
//     seen.insert(cur_vertex.clone());
    
//     while seen.len() < eord.len() {
//         let cur_edge = eord.get_mut(&cur_vertex).unwrap().remove(0);
//         cur_vertex = cur_edge[1..].to_string();
//         seen.insert(cur_vertex.clone());
//         t.push(cur_vertex.clone());
//     }
    
//     eord.get_mut(lst).unwrap().retain(|v| v != &dummy_vertex);
//     let mut new_seq = fst.to_string();
//     let mut cur_vertex = fst.to_string();
    
//     loop {
//         match eord.get_mut(&cur_vertex).and_then(|edges| edges.pop()) {
//             Some(cur_edge) => {
//                 cur_vertex = cur_edge[1..].to_string();
//                 new_seq.push(cur_edge.chars().last().unwrap());
//             }
//             None => {
//                 assert_eq!(cur_vertex, lst, "Last vertex should be equal to lst");
//                 break;
//             }
//         }
//     }
    
//     new_seq
// }

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
