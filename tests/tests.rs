use bioperm::altschul::klet_preserving_permutation;
use bioperm::utils::same_klets;
use quickcheck::quickcheck;

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
