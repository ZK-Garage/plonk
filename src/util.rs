// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use ark_ff::PrimeField;

/// Returns a vector of BlsScalars of increasing powers of x from x^0 to x^d.
pub(crate) fn powers_of<F: PrimeField>(
    scalar: &F,
    max_degree: usize,
) -> Vec<F> {
    let mut powers = Vec::with_capacity(max_degree + 1);
    powers.push(F::one());
    for i in 1..=max_degree {
        powers.push(powers[i - 1] * scalar);
    }
    powers
}
