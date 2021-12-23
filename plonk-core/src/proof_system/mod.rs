// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! PLONK Proving System

mod linearisation_poly;
mod permutation;
mod preprocess;
mod quotient_poly;
mod widget;

pub mod proof;
pub mod prover;
pub mod verifier;

pub use proof::*;
pub use prover::Prover;
pub use verifier::Verifier;
pub use widget::*;
