use super::poseidon_params::PoseidonParams;
use crate::merkle_tree::merkle_tree_fp::MerkleTreeHash;
use ark_ff::PrimeField;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct Poseidon<S: PrimeField> {
    pub(crate) params: Arc<PoseidonParams<S>>,
}

impl<S: PrimeField> Poseidon<S> {
    pub fn new(params: &Arc<PoseidonParams<S>>) -> Self {
        Poseidon {
            params: Arc::clone(params),
        }
    }

    pub fn get_t(&self) -> usize {
        self.params.t
    }

    pub fn permutation(&self, input: &[S]) -> Vec<S> {
        let t = self.params.t;
        assert_eq!(input.len(), t);

        let mut current_state = input.to_owned();
        for r in 0..self.params.rounds_f_beginning {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            current_state = self.matmul(&current_state, &self.params.mds);
        }
        let p_end = self.params.rounds_f_beginning + self.params.rounds_p;
        current_state = self.add_rc(&current_state, &self.params.opt_round_constants[0]);
        current_state = self.matmul(&current_state, &self.params.m_i);

        for r in self.params.rounds_f_beginning..p_end {
            current_state[0] = self.sbox_p(&current_state[0]);
            if r < p_end - 1 {
                current_state[0].add_assign(
                    &self.params.opt_round_constants[r + 1 - self.params.rounds_f_beginning][0],
                );
            }
            current_state = self.cheap_matmul(&current_state, p_end - r - 1);
        }
        for r in p_end..self.params.rounds {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            current_state = self.matmul(&current_state, &self.params.mds);
        }
        current_state
    }

    pub fn permutation_not_opt(&self, input: &[S]) -> Vec<S> {
        let t = self.params.t;
        assert_eq!(input.len(), t);

        let mut current_state = input.to_owned();

        for r in 0..self.params.rounds_f_beginning {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            current_state = self.matmul(&current_state, &self.params.mds);
        }
        let p_end = self.params.rounds_f_beginning + self.params.rounds_p;
        for r in self.params.rounds_f_beginning..p_end {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state[0] = self.sbox_p(&current_state[0]);
            current_state = self.matmul(&current_state, &self.params.mds);
        }
        for r in p_end..self.params.rounds {
            current_state = self.add_rc(&current_state, &self.params.round_constants[r]);
            current_state = self.sbox(&current_state);
            current_state = self.matmul(&current_state, &self.params.mds);
        }
        current_state
    }

    fn sbox(&self, input: &[S]) -> Vec<S> {
        input.iter().map(|el| self.sbox_p(el)).collect()
    }

    fn sbox_p(&self, input: &S) -> S {
        let mut input2 = *input;
        input2.square_in_place();

        match self.params.d {
            3 => {
                let mut out = input2;
                out.mul_assign(input);
                out
            }
            5 => {
                let mut out = input2;
                out.square_in_place();
                out.mul_assign(input);
                out
            }
            7 => {
                let mut out = input2;
                out.square_in_place();
                out.mul_assign(&input2);
                out.mul_assign(input);
                out
            }
            _ => {
                panic!()
            }
        }
    }

    fn cheap_matmul(&self, input: &[S], r: usize) -> Vec<S> {
        let v = &self.params.v[r];
        let w_hat = &self.params.w_hat[r];
        let t = self.params.t;

        let mut new_state = vec![S::zero(); t];
        new_state[0] = self.params.mds[0][0];
        new_state[0].mul_assign(&input[0]);
        for i in 1..t {
            let mut tmp = w_hat[i - 1];
            tmp.mul_assign(&input[i]);
            new_state[0].add_assign(&tmp);
        }
        for i in 1..t {
            new_state[i] = input[0];
            new_state[i].mul_assign(&v[i - 1]);
            new_state[i].add_assign(&input[i]);
        }

        new_state
    }

    fn matmul(&self, input: &[S], mat: &[Vec<S>]) -> Vec<S> {
        let t = mat.len();
        debug_assert!(t == input.len());
        let mut out = vec![S::zero(); t];
        for row in 0..t {
            for (col, inp) in input.iter().enumerate().take(t) {
                let mut tmp = mat[row][col];
                tmp.mul_assign(inp);
                out[row].add_assign(&tmp);
            }
        }
        out
    }

    fn add_rc(&self, input: &[S], rc: &[S]) -> Vec<S> {
        input
            .iter()
            .zip(rc.iter())
            .map(|(a, b)| {
                let mut r = *a;
                r.add_assign(b);
                r
            })
            .collect()
    }
}

impl<F: PrimeField> MerkleTreeHash<F> for Poseidon<F> {
    fn compress(&self, input: &[&F]) -> F {
        self.permutation(&[input[0].to_owned(), input[1].to_owned(), F::zero()])[0]
    }
}
#[cfg(test)]
mod poseidon_tests_bn256 {
    use super::*;
    use crate::fields::{bn256::FpBN256, utils::from_hex, utils::random_scalar};
    use crate::poseidon::poseidon_instance_bn256::POSEIDON_BN_PARAMS;

    type Scalar = FpBN256;

    static TESTRUNS: usize = 5;

    #[test]
    fn consistent_perm() {
        let poseidon = Poseidon::new(&POSEIDON_BN_PARAMS);
        let t = poseidon.params.t;
        for _ in 0..TESTRUNS {
            let input1: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

            let mut input2: Vec<Scalar>;
            loop {
                input2 = (0..t).map(|_| random_scalar()).collect();
                if input1 != input2 {
                    break;
                }
            }

            let perm1 = poseidon.permutation(&input1);
            let perm2 = poseidon.permutation(&input1);
            let perm3 = poseidon.permutation(&input2);
            assert_eq!(perm1, perm2);
            assert_ne!(perm1, perm3);
        }
    }

    #[test]
    fn kats() {
        let poseidon = Poseidon::new(&POSEIDON_BN_PARAMS);
        let input: Vec<Scalar> = vec![Scalar::from(0), Scalar::from(1), Scalar::from(2)];
        let perm = poseidon.permutation(&input);
        assert_eq!(
            perm[0],
            from_hex("0x2677d68d9cfa91f197bf5148b50afac461b6b8340ff119a5217794770baade5f")
        );
        assert_eq!(
            perm[1],
            from_hex("0x21ae9d716173496b62c76ad7deb4654961f64334441bcf77e17a047155a3239f")
        );
        assert_eq!(
            perm[2],
            from_hex("0x008f8e7c73ff20b6a141c48cef73215860acc749b14f0a7887f74950215169c6")
        );
    }
    #[test]
    fn opt_equals_not_opt() {
        let poseidon = Poseidon::new(&POSEIDON_BN_PARAMS);
        let t = poseidon.params.t;
        for _ in 0..TESTRUNS {
            let input: Vec<Scalar> = (0..t).map(|_| random_scalar()).collect();

            let perm1 = poseidon.permutation(&input);
            let perm2 = poseidon.permutation_not_opt(&input);
            assert_eq!(perm1, perm2);
        }
    }
}