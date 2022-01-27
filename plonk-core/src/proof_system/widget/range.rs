// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Range Gate

use crate::proof_system::{GateConstraint, GateValues};
use ark_ff::PrimeField;
use core::marker::PhantomData;

/// Range Gate
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct Range<F>(PhantomData<F>)
where
    F: PrimeField;

impl<F> GateConstraint<F> for Range<F>
where
    F: PrimeField,
{
    #[inline]
    fn constraints(separation_challenge: F, values: GateValues<F>) -> F {
        let four = F::from(4u64);
        let kappa = separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;
        let b_1 = delta(values.output - four * values.fourth);
        let b_2 = delta(values.right - four * values.output) * kappa;
        let b_3 = delta(values.left - four * values.right) * kappa_sq;
        let b_4 = delta(values.fourth_next - four * values.left) * kappa_cu;
        (b_1 + b_2 + b_3 + b_4) * separation_challenge
    }
}

/// Computes `f(f-1)(f-2)(f-3)`.
fn delta<F>(f: F) -> F
where
    F: PrimeField,
{
    let f_1 = f - F::one();
    let f_2 = f - F::from(2_u64);
    let f_3 = f - F::from(3_u64);
    f * f_1 * f_2 * f_3
}
