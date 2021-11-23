// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

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
