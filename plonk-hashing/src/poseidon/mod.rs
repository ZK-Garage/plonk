pub mod constants;
pub mod matrix;
pub mod mds;
pub mod poseidon_ref;
pub mod preprocessing;
pub mod round_constant;
pub mod round_numbers;
pub mod zprize_constraints;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum PoseidonError {
    #[error("Buffer is full")]
    FullBuffer,
}
