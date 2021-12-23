// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! Elliptic Curve Cryptography Gates

mod curve_addition;
mod fixed_base_scalar_mul;

pub use curve_addition::*;
pub use fixed_base_scalar_mul::*;
