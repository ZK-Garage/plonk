// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Collection of functions needed to use plonk library.
//!
//! Use this as the only import that you need to interact
//! with the principal data structures of the plonk library.

pub use crate::{
    circuit::{self, Circuit, PublicInputValue, VerifierData},
    constraint_system::{Point, StandardComposer, Variable},
    proof_system::{Prover, ProverKey, Verifier},
};

pub use crate::proof_system::{Proof, VerifierKey};

/// Collection of errors that the library exposes/uses.
pub use crate::error::Error;
