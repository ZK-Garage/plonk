// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! # PLONK
//! ![Build Status](https://github.com/rust-zkp/ark-plonk/workflows/Continuous%20integration/badge.svg)
//! [![Repository](https://img.shields.io/badge/github-plonk-blueviolet?logo=github)](https://github.com/rust-zkp/ark-plonk)
//! [![Documentation](https://img.shields.io/badge/docs-plonk-blue?logo=rust)](https://docs.rs/plonk/)
//!
//!
//! _This is a pure Rust implementation of the PLONK zk proving system_
//!
//!
//! ## About
//! Initial implementation created by [Kev](https://github.com/kevaundray), [Carlos](https://github.com/CPerezz) and [Luke](https://github.com/LukePearson1) at Dusk Network.
//! Redesigned by the [rust zkp](https://github.com/rust-zkp) team to have a backend which is compatible with the [arkworks](https://github.com/arkworks-rs) suite. This allows us to leverage the multitude of curves
//! and optimised algebra present in various arkworks repositories.
#![doc = include_str!("../test_circuit.md")]
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
pub mod prelude;
pub mod proof_system;

#[cfg(test)]
mod test;

#[doc = include_str!("../docs/notes-intro.md")]
pub mod notes {
    #[doc = include_str!("../docs/notes-commitments.md")]
    pub mod commitment_schemes {}
    #[doc = include_str!("../docs/notes-snark.md")]
    pub mod snark_construction {}
    #[doc = include_str!("../docs/notes-prove-verify.md")]
    pub mod prove_verify {}
    #[doc = include_str!("../docs/notes-KZG10.md")]
    pub mod kzg10_docs {}
}
