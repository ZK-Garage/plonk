// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Range Gate

use crate::proof_system::{
    linearisation_poly::CustomEvaluations, CustomValues, GateConstraint,
    WitnessValues,
};
use ark_ff::PrimeField;
use core::marker::PhantomData;

/// Values needed for the computation of the range gate constraint.
pub struct RangeVals<F>
where
    F: PrimeField,
{
    /// Fourth wire value in the next position
    pub d_next_val: F,
}

impl<F> CustomValues<F> for RangeVals<F>
where
    F: PrimeField,
{
    fn from_evaluations(custom_vals: &CustomEvaluations<F>) -> Self {
        let d_next_val = custom_vals.get("d_next_eval");
        RangeVals { d_next_val }
    }
}
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
    type CustomVals = RangeVals<F>;
    #[inline]
    fn constraints(
        separation_challenge: F,
        wit_vals: WitnessValues<F>,
        custom_vals: Self::CustomVals,
    ) -> F {
        let four = F::from(4u64);
        let kappa = separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;
        let b_1 = delta(wit_vals.c_val - four * wit_vals.d_val);
        let b_2 = delta(wit_vals.b_val - four * wit_vals.c_val) * kappa;
        let b_3 = delta(wit_vals.a_val - four * wit_vals.b_val) * kappa_sq;
        let b_4 =
            delta(custom_vals.d_next_val - four * wit_vals.a_val) * kappa_cu;
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
