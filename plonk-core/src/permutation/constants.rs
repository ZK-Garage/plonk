// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! Constants used in the permutation argument to ensure that the wire subsets
//! are disjoint.

#![allow(non_snake_case)]

use ark_ff::PrimeField;

pub(crate) fn K1<F: PrimeField>() -> F {
    F::from(7_u64)
}

pub(crate) fn K2<F: PrimeField>() -> F {
    F::from(13_u64)
}

pub(crate) fn K3<F: PrimeField>() -> F {
    F::from(17_u64)
}
