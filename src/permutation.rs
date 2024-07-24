use crate::spec::{Spec, State};
use halo2_proofs::arithmetic::FieldExt;

impl<const T: usize, const RATE: usize> Spec<T, RATE> {
    /// Applies the Poseidon2 permutation to the given state
    pub fn permute(&self, state: &mut State<T>) {
        let (r_f, r_p) = (self.r_f / 2, self.r_p);

        self.mds_external.apply(state);
        for constants in self.constants.iter().take(r_f) {
            state.add_constants(constants);
            state.sbox_full();
            self.mds_external.apply(state);
        }

        for constants in self.constants.iter().skip(r_f).take(r_p) {
            state.add_constants(constants);
            state.sbox_part();
            self.mds_internal.apply(state);
        }

        for constants in self.constants.iter().skip(r_f + r_p) {
            state.add_constants(constants);
            state.sbox_full();
            self.mds_external.apply(state);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::State;
    use crate::poseidon2_instance::from_hex;
    use crate::spec::Spec;
    use halo2_proofs::pairing::bn256::Fr;
    use halo2_proofs::pairing::group::ff::PrimeField;

    #[test]
    fn poseidon2_test() {
        use halo2_proofs::pairing::group::ff::Field;
        use rand_core::OsRng;
        use std::time::Instant;

        macro_rules! run_test {
            (
                $([$RF:expr, $RP:expr, $T:expr, $RATE:expr]),*
            ) => {
                $(
                    {
                        const R_F: usize = $RF;
                        const R_P: usize = $RP;
                        const T: usize = $T;
                        const RATE: usize = $RATE;
                        let mut state = State(
                            (0..T)
                                .map(|_| Fr::random(OsRng))
                                .collect::<Vec<Fr>>()
                                .try_into().unwrap(),
                        );

                        let spec = Spec::<T, RATE>::new(R_F, R_P);
                        let now = Instant::now();
                        {
                            spec.permute(&mut state);
                        }
                        let elapsed = now.elapsed();
                        println!("Elapsed: {:.2?}", elapsed);
                    }
                )*
            };
        }
        run_test!([8, 56, 3, 2]);
        // run_test!([8, 57, 4, 3]);
        // run_test!([8, 57, 5, 4]);
        // run_test!([8, 57, 6, 5]);
        // run_test!([8, 57, 7, 6]);
        // run_test!([8, 57, 8, 7]);
        // run_test!([8, 57, 9, 8]);
        // run_test!([8, 57, 10, 9]);
    }
    #[test]
    fn test_spec() {
        const R_F: usize = 8;
        const R_P: usize = 56;
        const T: usize = 3;
        const RATE: usize = 2;

        let state: State<3> = State(
            vec![0u64, 1, 2]
                .into_iter()
                .map(Fr::from)
                .collect::<Vec<Fr>>()
                .try_into()
                .unwrap(),
        );
        assert_eq!(state.words()[0], from_hex("0"));
        assert_eq!(state.words()[1], from_hex("1"));
        assert_eq!(state.words()[2], from_hex("2"));

        let spec = Spec::<T, RATE>::new_from_instance(R_F, R_P);

        assert_eq!(spec.constants.len(), R_F + R_P);
        assert_eq!(
            spec.constants[0],
            [
                from_hex("0x1d066a255517b7fd8bddd3a93f7804ef7f8fcde48bb4c37a59a09a1a97052816"),
                from_hex("0x29daefb55f6f2dc6ac3f089cebcc6120b7c6fef31367b68eb7238547d32c1610"),
                from_hex("0x1f2cb1624a78ee001ecbd88ad959d7012572d76f08ec5c4f9e8b7ad7b0b4e1d1")
            ]
        );
        assert_eq!(
            spec.constants[4],
            [
                from_hex("0x1a1d063e54b1e764b63e1855bff015b8cedd192f47308731499573f23597d4b5"),
                from_hex("0x0000000000000000000000000000000000000000000000000000000000000000"),
                from_hex("0x0000000000000000000000000000000000000000000000000000000000000000"),
            ]
        );

        assert_eq!(
            spec.mds_external.0 .0[0],
            [
                Fr::from_str_vartime("2").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
            ]
        );
        assert_eq!(
            spec.mds_external.0 .0[1],
            [
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("2").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
            ]
        );
        assert_eq!(
            spec.mds_external.0 .0[2],
            [
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("2").unwrap(),
            ]
        );

        assert_eq!(
            spec.mds_internal.0 .0[0],
            [
                Fr::from_str_vartime("2").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
            ]
        );
        assert_eq!(
            spec.mds_internal.0 .0[1],
            [
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("2").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
            ]
        );
        assert_eq!(
            spec.mds_internal.0 .0[2],
            [
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("3").unwrap(),
            ]
        );
    }

    #[test]
    fn test_against_test_vectors() {
        // https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/test_vectors.txt
        // poseidon2perm_x5_254_3
        const R_F: usize = 8;
        const R_P: usize = 56;
        const T: usize = 3;
        const RATE: usize = 2;

        let state: State<3> = State(
            vec![0u64, 1, 2]
                .into_iter()
                .map(Fr::from)
                .collect::<Vec<Fr>>()
                .try_into()
                .unwrap(),
        );
        assert_eq!(state.words()[0], from_hex("0"));
        assert_eq!(state.words()[1], from_hex("1"));
        assert_eq!(state.words()[2], from_hex("2"));

        let spec = Spec::<T, RATE>::new_from_instance(R_F, R_P);

        let mut state_0 = state;
        spec.permute(&mut state_0);
        let expected = vec![
            "0bb61d24daca55eebcb1929a82650f328134334da98ea4f847f760054f4a3033",
            "303b6f7c86d043bfcbcc80214f26a30277a15d3f74ca654992defe7ff8d03570",
            "1ed25194542b12eef8617361c3ba7c52e660b145994427cc86296242cf766ec8",
        ];

        for (word, expected) in state_0.words().into_iter().zip(expected.iter()) {
            assert_eq!(word, from_hex::<Fr>(expected));
        }
    }
}
