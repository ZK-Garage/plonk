// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Collection of functions needed to build and use PLONK circuits.
//!
//! Use this as the only import that you need to interact with the principal
//! data structures of the plonk library.

pub use crate::{
    circuit::{self, verify_proof, Circuit, VerifierData},
    constraint_system::{ecc::Point, StandardComposer, Variable},
    error::Error,
    proof_system::{Proof, ProverKey, VerifierKey},
    util::from_embedded_curve_scalar,
};
