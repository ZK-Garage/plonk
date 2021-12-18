// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE 
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-INFRA. All rights reserved.

//! The constraint system module stores the implementation of the PLONK
//! [`StandardComposer`], as well as the circuit tools and abstractions, used by
//! the Composer to generate, build, preprocess circuits.

mod arithmetic;
mod boolean;
mod logic;
mod range;

pub(crate) mod composer;
pub(crate) mod helper;
pub(crate) mod variable;

pub mod ecc;

pub(crate) use variable::WireData;

pub use composer::StandardComposer;
pub use variable::Variable;
