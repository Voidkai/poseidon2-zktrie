
#[allow(clippy::too_many_arguments)]
#[allow(clippy::derive_hash_xor_eq)]
pub mod bn256;

pub mod utils;

// sage:
// p = 21888242871839275222246405745257275088548364400416034343698204186575808495617
// F = GF(p)
// F.multiplicative_generator()
// F(7).is_primitive_root()
