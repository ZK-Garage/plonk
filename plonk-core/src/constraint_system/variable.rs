// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! This module holds the components needed in the Constraint System.
//!
//! The two components used are Variables and Wires.
use std::fmt::Display;

/// The value is a reference to the actual value that was added to the
/// constraint system
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Variable(pub(crate) usize);

impl Display for Variable {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Stores the data for a specific wire in an arithmetic circuit
/// This data is the gate index and the type of wire
/// Left(1) signifies that this wire belongs to the first gate and is the left
/// wire
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum WireData {
    /// Left Wire of n'th gate
    Left(usize),
    /// Right Wire of n'th gate
    Right(usize),
    /// Output Wire of n'th gate
    Output(usize),
    /// Fourth Wire of n'th gate
    Fourth(usize),
}
