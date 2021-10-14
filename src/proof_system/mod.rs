// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Proving system
mod aggregated_proof;
mod preprocess;
/// Represents a PLONK Prover
pub mod prover;
pub use proof::*;
pub(crate) mod quotient_poly;
/// Represents a PLONK Verifier
pub mod verifier;
pub use prover::Prover;
pub use verifier::Verifier;
pub use widget::*;

pub(crate) use aggregated_proof::PCAggregateProof;
/// Represents a PLONK Proof
pub mod proof;
pub(crate) mod widget;
pub use proof::Proof;
pub use widget::VerifierKey;
pub(crate) mod linearisation_poly;
