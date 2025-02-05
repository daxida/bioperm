use std::collections::{HashMap, HashSet};

use bioperm::altschul::klet_preserving_permutation;
use bioperm::utils::same_klets;
use quickcheck::quickcheck;

fn klet_uniform<F>(seq: &str, all_perms: &[&str], k: usize, tested_fn: F)
where
    F: Fn(&str, usize, bool) -> String,
{
    let all_perms_set: HashSet<_> = all_perms.iter().map(|s| s.to_string()).collect();

    let mut cnt: HashMap<String, usize> = HashMap::new();
    let max_iterations = 1000;
    for _ in 0..max_iterations {
        let res = tested_fn(seq, k, false);
        *cnt.entry(res.clone()).or_insert(0) += 1;
        assert!(same_klets(seq, &res, k));
    }

    for (perm, &amount) in &cnt {
        assert!(all_perms_set.contains(perm));
        let ratio = amount as f64 / max_iterations as f64 - 1.0 / all_perms_set.len() as f64;
        assert!(ratio.abs() < 0.05);
    }

    let remaining: HashSet<_> = cnt.keys().cloned().collect();
    assert!(remaining.is_subset(&all_perms_set));
}

#[test]
fn test_klet_uniform_altschul() {
    let seq = "ACTAGTAT";
    let all_perms = [
        "AGTACTAT", "ACTAGTAT", "ATAGTACT", "ATACTAGT", "AGTATACT", "ACTATAGT",
    ];
    klet_uniform(seq, &all_perms, 2, klet_preserving_permutation);

    let seq = "AAATAAA";
    let all_perms = ["AAAATAA", "AATAAAA", "AAATAAA"];
    klet_uniform(seq, &all_perms, 3, klet_preserving_permutation);

    let seq = "AAAATAAAA";
    let all_perms = ["AAAAATAAA", "AAATAAAAA", "AAAATAAAA"];
    klet_uniform(seq, &all_perms, 4, klet_preserving_permutation);
}

// Includes the k in order to only make tests such that k < len(seq)
#[derive(Debug, Clone)]
struct SequenceK {
    seq: String,
    k: usize,
}

impl quickcheck::Arbitrary for SequenceK {
    fn arbitrary(g: &mut quickcheck::Gen) -> Self {
        let k = usize::arbitrary(g) % 5 + 2;
        let wlen = usize::arbitrary(g) % 200 + 1 + k;
        let letters = ['A', 'C', 'G', 'T'];
        let seq: String = (0..wlen).map(|_| *g.choose(&letters).unwrap()).collect();
        Self { seq, k }
    }
}

quickcheck! {
    fn test_klet_preserving_permutation(seqk: SequenceK) -> bool {
        let res = klet_preserving_permutation(&seqk.seq, seqk.k, false);
        (2..=seqk.k).all(|j| same_klets(&seqk.seq, &res, j))
    }
}
