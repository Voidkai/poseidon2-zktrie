//! Poseidon hashing implemention with variable lenght input setting. This crate
//! also exposes constant parameters for circuit implementations

#![deny(missing_debug_implementations)]
#![deny(missing_docs)]

//mod grain;
mod matrix;
mod permutation;
mod poseidon2;
mod poseidon2_instance;
mod spec;

pub use crate::poseidon2::Poseidon2;
pub use crate::spec::{MDSMatrix, SparseMDSMatrix, Spec, State};
