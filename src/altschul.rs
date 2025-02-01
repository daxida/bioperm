use crate::unionfind::UnionFind;
use crate::utils::{each_cons, same_klets};
use rand::seq::SliceRandom;
use std::collections::{HashMap, HashSet};

pub const S1: &str = "AGACATAAAGTTCCGTACTGCCGGGAT";

fn check_if_connected(z_graph: &[&str], vertices: &[&str]) -> bool {
    debug_assert_eq!(z_graph.len(), vertices.len() - 1);

    let mut uf = UnionFind::new(vertices.len());
    for edge in z_graph {
        let fr = &edge[..edge.len() - 1];
        let to = &edge[1..];
        if let (Some(fr_idx), Some(to_idx)) = (
            vertices.iter().position(|v| *v == fr),
            vertices.iter().position(|v| *v == to),
        ) {
            uf.union(fr_idx, to_idx);
        } else {
            panic!();
        }
    }

    uf.group_count() == 1
}

fn edge_ordering(seq: &str, k: usize) -> HashMap<&str, Vec<&str>> {
    let mut eord: HashMap<&str, Vec<&str>> = HashMap::new();

    for chunk in each_cons(seq, k) {
        let key = &chunk[..k - 1];
        eord.entry(key).or_default().push(chunk);
    }

    // Ensure the last vertex is added
    let last_vertex = &seq[seq.len() - (k - 1)..];
    eord.entry(last_vertex).or_default();

    eord
}

pub fn klet_preserving_permutation(seq: &str, k: usize, debug: bool) -> String {
    assert!(1 < k);
    assert!(k < seq.len());

    let mut eord = edge_ordering(seq, k);
    let s1 = &seq[..k - 1];
    let sf = &seq[seq.len() - (k - 1)..];

    if debug {
        println!("1. {:?}", seq);
        println!("Edge ordering: {:?}", eord);
        println!("s1: {:?}, sf: {:?}", s1, sf);
    }

    let mut last_edges: Vec<&str> = Vec::new();
    let vertices: Vec<&str> = eord.keys().cloned().collect();
    let mut seen = HashSet::new();
    let mut max_iterations = 1000;

    while max_iterations > 0 {
        last_edges.clear();
        for (vertex, edges) in &eord {
            if *vertex != sf {
                if let Some(edge) = edges.choose(&mut rand::thread_rng()) {
                    last_edges.push(edge);
                }
            }
        }

        let key = last_edges.clone();
        if !seen.insert(key) {
            continue;
        }
        if check_if_connected(&last_edges, &vertices) {
            break;
        }
        max_iterations -= 1;
    }

    if max_iterations == 0 {
        panic!("Exhausted iterations. Returning original sequence.");
    }

    if debug {
        println!("5. {:?}", last_edges);
    }

    for last_edge in &last_edges {
        let vertex = &last_edge[..k - 1];
        if let Some(edges) = eord.get_mut(vertex) {
            if let Some(pos) = edges.iter().position(|e| e == last_edge) {
                edges.remove(pos);
                edges.push(last_edge);
            }
        }
    }

    let mut eord_perm = HashMap::new();
    for (vertex, mut edges) in eord {
        if vertex == sf {
            edges.shuffle(&mut rand::thread_rng());
            eord_perm.insert(vertex, edges);
        } else {
            let mut to_permute = edges[..edges.len() - 1].to_vec();
            to_permute.shuffle(&mut rand::thread_rng());
            to_permute.push(edges[edges.len() - 1]);
            eord_perm.insert(vertex, to_permute);
        }
    }

    let mut new_seq = String::from(s1);
    let mut cur_vertex = s1.to_string();
    loop {
        match eord_perm.get_mut(&cur_vertex[..]) {
            Some(edges) if !edges.is_empty() => {
                let cur_edge = edges.remove(0);
                cur_vertex = cur_edge[1..].to_string();
                new_seq.push(cur_edge.chars().last().unwrap());
                if debug {
                    println!("{:?}, {:?}, {:?}", cur_edge, cur_vertex, new_seq);
                }
            }
            _ => {
                assert_eq!(cur_vertex, sf);
                break;
            }
        }
    }

    for (key, value) in &eord_perm {
        if !value.is_empty() {
            panic!("{} edge list was not exhausted", key);
        }
    }

    assert_eq!(seq.len(), new_seq.len(), "Sequence length mismatch");
    assert!(same_klets(seq, &new_seq, 2), "Doublets are not preserved");

    new_seq
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_if_connected() {
        let z_graph = ["AT", "CT", "GA"];
        let vertices = ["A", "C", "G", "T"];
        assert!(check_if_connected(&z_graph, &vertices));
    }

    #[test]
    fn test_edge_ordering() {
        let result = edge_ordering(S1, 3);
        assert_eq!(result.len(), 16);
    }

    #[test]
    fn test_klet_preservation() {
        let _seq = "AGTACTAT".repeat(2);
        let seq = _seq.as_str();
        let seq_2 = klet_preserving_permutation(seq, 2, false);
        assert!(same_klets(seq, &seq_2, 2));

        let seq_3 = klet_preserving_permutation(seq, 3, false);
        assert!(same_klets(seq, &seq_3, 2));
        assert!(same_klets(seq, &seq_3, 3));
    }
}
