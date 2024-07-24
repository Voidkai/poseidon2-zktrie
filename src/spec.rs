use std::ops::AddAssign;
use std::ops::Index;
use std::ops::MulAssign;
//use crate::grain::Grain;
use crate::matrix::Matrix;
use crate::poseidon2_instance::RC3;
use halo2_proofs::arithmetic::Field;
use halo2_proofs::pairing::group::ff::PrimeField;
use halo2_proofs::{arithmetic::FieldExt, pairing::bn256::Fr};

/// `State` is structure `T` sized field elements that are subjected to
/// permutation
#[derive(Clone, Debug, PartialEq)]
pub struct State<const T: usize>(pub(crate) [Fr; T]);

impl<const T: usize> Default for State<T> {
    /// The capacity value is 2**64 + (o − 1) where o the output length.
    fn default() -> Self {
        let mut state = [Fr::zero(); T];
        state[0] = Fr::from_u128(1 << 64);
        State(state)
    }
}

impl<const T: usize> State<T> {
    /// Applies sbox for all elements of the state.
    /// Only supports `alpha = 5` sbox case.
    pub(crate) fn sbox_full(&mut self) {
        for e in self.0.iter_mut() {
            let tmp = e.mul(e);
            e.mul_assign(tmp);
            e.mul_assign(tmp);
        }
    }

    /// Partial round sbox applies sbox to the first element of the state.
    /// Only supports `alpha = 5` sbox case
    pub(crate) fn sbox_part(&mut self) {
        let tmp = self.0[0].mul(&self.0[0]);
        self.0[0].mul_assign(tmp);
        self.0[0].mul_assign(tmp);
    }

    /// Adds constants to all elements of the state
    pub(crate) fn add_constants(&mut self, constants: &[Fr; T]) {
        for (e, constant) in self.0.iter_mut().zip(constants.iter()) {
            e.add_assign(constant)
        }
    }

    /// Only adds a constant to the first element of the state.It is used with
    /// optimized rounds constants where only single element is added in
    /// each partial round
    pub(crate) fn add_constant(&mut self, constant: &Fr) {
        self.0[0].add_assign(constant)
    }

    /// Copies elements of the state
    pub fn words(&self) -> [Fr; T] {
        self.0
    }

    /// Second element of the state is the result
    pub(crate) fn result(&self) -> Fr {
        self.0[1]
    }
}

/// `Spec` holds construction parameters as well as constants that are used in
/// permutation step. Constants are planned to be hardcoded once transcript
/// design matures. Number of partial rounds can be deriven from number of
/// constants.
#[derive(Debug, Clone)]
pub struct Spec<const T: usize, const RATE: usize> {
    pub(crate) r_f: usize,
    pub(crate) r_p: usize,
    pub(crate) mds_external: MDSMatrix<T, RATE>,
    pub(crate) mds_internal: MDSMatrix<T, RATE>,
    pub(crate) constants: Vec<[Fr; T]>,
}

// impl<F: FieldExt, const T: usize, const RATE: usize> Spec<F, T, RATE> {
//     /// Number of full rounds
//     pub fn r_f(&self) -> usize {
//         self.r_f.clone()
//     }
//     /// Set of MDS Matrices used in permutation line
//     pub fn mds_matrices(&self) -> &MDSMatrices<F, T, RATE> {
//         &self.mds_matrices
//     }
//     /// Optimised round constants
//     pub fn constants(&self) -> &OptimizedConstants<F, T> {
//         &self.constants
//     }
// }

/// `MDSMatrix` is applied to `State` to achive linear layer of Poseidon
#[derive(Clone, Debug)]
pub struct MDSMatrix<const T: usize, const RATE: usize>(pub(crate) Matrix<Fr, T>);

impl<const T: usize, const RATE: usize> Index<usize> for MDSMatrix<T, RATE> {
    type Output = [Fr; T];

    fn index(&self, idx: usize) -> &Self::Output {
        &self.0 .0[idx]
    }
}

impl<const T: usize, const RATE: usize> MDSMatrix<T, RATE> {
    /// Applies `MDSMatrix` to the state
    pub(crate) fn apply(&self, state: &mut State<T>) {
        state.0 = self.0.mul_vector(&state.0);
    }

    /// Given two `T` sized vector constructs the `t * t` Cauchy matrix
    pub(super) fn cauchy(xs: &[Fr; T], ys: &[Fr; T]) -> Self {
        let mut m = Matrix::default();
        for (i, x) in xs.iter().enumerate() {
            for (j, y) in ys.iter().enumerate() {
                let sum = *x + *y;
                debug_assert!(!sum.is_zero_vartime());
                m.set(i, j, sum.invert().unwrap());
            }
        }
        MDSMatrix(m)
    }

    /// Inverts the MDS matrix
    fn invert(&self) -> Self {
        Self(self.0.invert())
    }

    /// Used in calculation of optimized round constants. Calculates `v' = M *
    /// v` where vectors are `T` sized
    fn mul_constants(&self, v: &[Fr; T]) -> [Fr; T] {
        self.0.mul_vector(v)
    }

    /// Multiplies two MDS matrices. Used in sparse matrix calculations
    fn mul(&self, other: &Self) -> Self {
        Self(self.0.mul(&other.0))
    }

    fn transpose(&self) -> Self {
        Self(self.0.transpose())
    }

    /// See Section B in Supplementary Material https://eprint.iacr.org/2019/458.pdf
    /// Factorises an MDS matrix `M` into `M'` and `M''` where `M = M' *  M''`.
    /// Resulted `M''` matrices are the sparse ones while `M'` will contribute
    /// to the accumulator of the process
    fn factorise(&self) -> (Self, SparseMDSMatrix<T, RATE>) {
        // Given `(t-1 * t-1)` MDS matrix called `hat` constructs the matrix in
        // form `[[1 | 0], [0 | m]]`
        let prime = |hat: Matrix<Fr, RATE>| -> MDSMatrix<T, RATE> {
            let mut prime = Matrix::identity();
            for (prime_row, hat_row) in prime.0.iter_mut().skip(1).zip(hat.0.iter()) {
                for (el_prime, el_hat) in prime_row.iter_mut().skip(1).zip(hat_row.iter()) {
                    *el_prime = *el_hat;
                }
            }
            Self(prime)
        };

        // Given `(t-1)` sized `w_hat` vector constructs the matrix in form
        // `[[m_0_0 | m_0_i], [w_hat | identity]]`
        let prime_prime = |w_hat: [Fr; RATE]| -> Self {
            let mut prime_prime = Matrix::identity();
            prime_prime.0[0] = self.0 .0[0];
            for (row, w) in prime_prime.0.iter_mut().skip(1).zip(w_hat.iter()) {
                row[0] = *w
            }
            Self(prime_prime)
        };

        let w = self.0.w();
        let m_hat = self.0.sub::<RATE>();
        let m_hat_inverse = m_hat.invert();
        let w_hat = m_hat_inverse.mul_vector(&w);
        (prime(m_hat), prime_prime(w_hat).transpose().into())
    }

    /// Returns rows of the MDS matrix
    pub fn rows(&self) -> [[Fr; T]; T] {
        self.0 .0
    }
}

/// `SparseMDSMatrix` are in `[row], [hat | identity]` form and used in linear
/// layer of partial rounds instead of the original MDS
#[derive(Debug, Clone)]
pub struct SparseMDSMatrix<const T: usize, const RATE: usize> {
    pub(crate) row: [Fr; T],
    pub(crate) col_hat: [Fr; RATE],
}

impl<const T: usize, const RATE: usize> SparseMDSMatrix<T, RATE> {
    /// Returns the first row
    pub fn row(&self) -> &[Fr; T] {
        &self.row
    }

    /// Returns the first column without first element in the first row
    pub fn col_hat(&self) -> &[Fr; RATE] {
        &self.col_hat
    }

    /// Applies the sparse MDS matrix to the state
    pub(crate) fn apply(&self, state: &mut State<T>) {
        let words = state.words();
        state.0[0] = self
            .row
            .iter()
            .zip(words.iter())
            .fold(Fr::zero(), |acc, (e, cell)| acc + (*e * *cell));

        for ((new_word, col_el), word) in (state.0)
            .iter_mut()
            .skip(1)
            .zip(self.col_hat.iter())
            .zip(words.iter().skip(1))
        {
            *new_word = *col_el * words[0] + word;
        }
    }
}

impl<const T: usize, const RATE: usize> From<MDSMatrix<T, RATE>> for SparseMDSMatrix<T, RATE> {
    /// Assert the form and represent an MDS matrix as a sparse MDS matrix
    fn from(mds: MDSMatrix<T, RATE>) -> Self {
        let mds = mds.0;
        for (i, row) in mds.0.iter().enumerate().skip(1) {
            for (j, _) in row.iter().enumerate().skip(1) {
                assert_eq!(row[j], if i != j { Fr::zero() } else { Fr::one() });
            }
        }

        let (mut row, mut col_hat) = ([Fr::zero(); T], [Fr::zero(); RATE]);
        for (row_el, el) in row.iter_mut().zip(mds.0[0].iter()) {
            *row_el = *el
        }
        for (col_el, row) in col_hat.iter_mut().zip(mds.0.iter().skip(1)) {
            *col_el = row[0]
        }

        SparseMDSMatrix { row, col_hat }
    }
}

impl<const T: usize, const RATE: usize> Spec<T, RATE> {
    /// Given number of round parameters constructs new Posedion instance
    /// calculating unoptimized round constants with reference `Grain` then
    /// calculates optimized constants and sparse matrices
    pub fn new(r_f: usize, r_p: usize) -> Self {
        let constants = (0..(r_f + r_p))
            .map(|_| {
                let mut round_constants = [Fr::zero(); T];
                for c in round_constants.iter_mut() {
                    *c = Fr::zero();
                }
                round_constants
            })
            .collect::<Vec<[Fr; T]>>();
        let mds = MDSMatrix::<T, RATE>(Matrix::identity());
        let mds_external = mds.clone();
        let mds_internal = mds.clone();
        Self {
            r_f,
            r_p,
            mds_external,
            mds_internal,
            constants,
        }
    }
}

impl Spec<3, 2> {
    /// Given number of round parameters constructs new Posedion instance
    /// calculating unoptimized round constants with reference `Grain` then
    /// calculates optimized constants and sparse matrices
    pub fn new_from_instance(r_f: usize, r_p: usize) -> Self {
        // let (constants, mds) = Grain::generate(r_f, r_p);
        let constants = (0..(r_f + r_p))
            .map(|i| {
                let mut round_constants = [Fr::zero(); 3];
                for (j, c) in round_constants.iter_mut().enumerate() {
                    *c = RC3[i][j];
                }
                round_constants
            })
            .collect::<Vec<[Fr; 3]>>();

        let mds_external = MDSMatrix::<3, 2>(Matrix([
            [
                Fr::from_str_vartime("2").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
            ],
            [
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("2").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
            ],
            [
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("2").unwrap(),
            ],
        ]));

        let mds_internal = MDSMatrix::<3, 2>(Matrix([
            [
                Fr::from_str_vartime("2").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
            ],
            [
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("2").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
            ],
            [
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("1").unwrap(),
                Fr::from_str_vartime("3").unwrap(),
            ],
        ]));

        Self {
            r_f,
            r_p,
            mds_external,
            mds_internal,
            constants,
        }
    }
}
