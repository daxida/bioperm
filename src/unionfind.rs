pub struct UnionFind {
    parents: Vec<i32>,
    group: usize,
}

impl UnionFind {
    pub fn new(n: usize) -> UnionFind {
        UnionFind {
            parents: vec![-1; n],
            group: n,
        }
    }

    pub fn find(&mut self, x: usize) -> i32 {
        if self.parents[x] < 0 {
            x as i32
        } else {
            self.parents[x] = self.find(self.parents[x] as usize);
            self.parents[x]
        }
    }

    pub fn union(&mut self, x: usize, y: usize) {
        let mut new_x = self.find(x);
        let mut new_y = self.find(y);

        if new_x != new_y {
            self.group -= 1;
            if self.parents[new_x as usize] > self.parents[new_y as usize] {
                std::mem::swap(&mut new_y, &mut new_x);
            }
            self.parents[new_x as usize] += self.parents[new_y as usize];
            self.parents[new_y as usize] = new_x;
        }
    }

    pub fn group_count(&self) -> usize {
        self.group
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basics() {
        let mut uf = UnionFind::new(3);
        uf.union(1, 1);
        assert!(uf.group_count() == 3);
        uf.union(1, 2);
        assert!(uf.group_count() == 2);
        uf.union(0, 2);
        assert!(uf.group_count() == 1);
    }
}
