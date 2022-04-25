// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Permutations over Lagrange-bases for Oecumenical Noninteractive
//! arguments of Knowledge (PLONK) is a zero knowledge proof system.
//!
//! This protocol was created by:
//! - Ariel Gabizon (Protocol Labs),
//! - Zachary J. Williamson (Aztec Protocol)
//! - Oana Ciobotaru
//!
//! This crate contains a pure Rust implementation of this algorithm using
//! code done by the creators of the protocol as a reference implementation:
//!
//! <https://github.com/AztecProtocol/barretenberg/blob/master/barretenberg/src/aztec/plonk/>

// Bitshift/Bitwise ops are allowed to gain performance.
#![allow(clippy::suspicious_arithmetic_impl)]
// Some structs do not have AddAssign or MulAssign impl.
#![allow(clippy::suspicious_op_assign_impl)]
// Variables have always the same names in respect to wires.
#![allow(clippy::many_single_char_names)]
// Bool expr are usually easier to read with match statements.
#![allow(clippy::match_bool)]
// We have quite some functions that require quite some args by it's nature.
// It can be refactored but for now, we avoid these warns.
#![allow(clippy::too_many_arguments)]
#![deny(rustdoc::broken_intra_doc_links)]
#![deny(missing_docs)]

extern crate alloc;

mod permutation;
mod transcript;
mod util;

pub mod circuit;
pub mod commitment;
pub mod constraint_system;
pub mod error;
pub mod lookup;
pub mod prelude;
pub mod proof_system;

#[cfg(test)]
mod test;
