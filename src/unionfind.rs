pub struct UnionFind {
    n: i32,
    parents: Vec<i32>,
    group: usize,
}

impl UnionFind {
    fn new(n: i32) -> UnionFind {
        UnionFind {
            n: n,
            parents: vec![-1; n as usize],
            group: n as usize,
        }
    }

    fn find(&mut self, x: i32) -> i32 {
        if self.parents[x as usize] < 0 {
            return x;
        } else {
            self.parents[x as usize] = self.find(self.parents[x as usize]);
            return self.parents[x as usize];
        }
    }

    fn union(&mut self, x: i32, y: i32) {
        let mut new_x = self.find(x);
        let mut new_y = self.find(y);

        if new_x != new_y {
            self.group -= 1;
            if self.parents[new_x as usize] > self.parents[new_y as usize] {
                let tmp: i32 = new_y;
                new_y = new_x;
                new_x = tmp;
            }
            self.parents[new_x as usize] += self.parents[new_y as usize];
            self.parents[new_y as usize] = new_x;
        }
    }

    fn group_count(&self) -> usize {
        return self.group;
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
