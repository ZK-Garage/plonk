// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! This module holds the components needed in the Constraint System.
//!
//! The two components used are Variables and Wires.

/// The value is a reference to the actual value that was added to the
/// constraint system
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub struct Variable(pub(crate) usize);

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
