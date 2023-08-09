// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! The constraint system module stores the implementation of the PLONK
//! [`StandardComposer`], as well as the circuit tools and abstractions, used by
//! the Composer to generate, build, preprocess circuits.

mod arithmetic;
mod boolean;
mod hash;
mod logic;
mod lookup;
mod range;

pub(crate) mod composer;
pub(crate) mod helper;
pub(crate) mod variable;

pub mod ecc;

pub(crate) use hash::SBOX_ALPHA;
pub(crate) use variable::WireData;

pub use composer::StandardComposer;
pub use variable::Variable;
