#![allow(non_snake_case)]
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use ark_ff::PrimeField;

/// Constants used in the permutation argument to ensure that the wire subsets
/// are disjoint.
pub(crate) fn K1<F: PrimeField>() -> F {
    F::from(7 as u64)
}
pub(crate) fn K2<F: PrimeField>() -> F {
    F::from(13 as u64)
}
pub(crate) fn K3<F: PrimeField>() -> F {
    F::from(17 as u64)
}
