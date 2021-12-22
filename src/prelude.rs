// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-INFRA. All rights reserved.

//! Collection of functions needed to build and use PLONK circuits.
//!
//! Use this as the only import that you need to interact with the principal
//! data structures of the plonk library.

pub use crate::{
    circuit::{
        self, verify_proof, Circuit, FeIntoPubInput, GeIntoPubInput,
        PublicInputValue, VerifierData,
    },
    constraint_system::{ecc::Point, StandardComposer, Variable},
    error::Error,
    proof_system::{Proof, ProverKey, VerifierKey},
    util::from_embedded_curve_scalar,
};
