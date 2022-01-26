// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved..

//! Logic Gates

<<<<<<< HEAD
use crate::proof_system::widget::{GateConstraint, WitnessValues};
use ark_ff::PrimeField;
=======
use crate::{
    get_label,
    proof_system::{
        linearisation_poly::CustomEvaluations,
        widget::{GateConstraint, WitnessValues},
    },
};
use ark_ff::Field;
>>>>>>> 8d683a0 (Add `from_evaluations` constructor)
use core::marker::PhantomData;

use super::CustomValues;

pub struct LogicVals<F>
where
    F: PrimeField,
{
    pub a_next_eval: F,
    pub b_next_eval: F,
    pub d_next_eval: F,
    pub q_c_eval: F,
}

impl<F> CustomValues<F> for LogicVals<F>
where
    F: Field,
{
    fn from_evaluations(custom_evals: CustomEvaluations<F>) -> Self {
        let a_next_eval = custom_evals.get(get_label!(a_next_eval));
        let b_next_eval = custom_evals.get(get_label!(b_next_eval));
        let d_next_eval = custom_evals.get(get_label!(d_next_eval));
        let q_c_eval = custom_evals.get(get_label!(q_c_eval));
        LogicVals {
            a_next_eval,
            b_next_eval,
            d_next_eval,
            q_c_eval,
        }
    }
}

<<<<<<< HEAD
impl<F> CustomValues<F> for LogicVals<F> where F: PrimeField {}
=======
>>>>>>> 8d683a0 (Add `from_evaluations` constructor)
/// Logic Gate
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct Logic<F>(PhantomData<F>)
where
    F: PrimeField;

impl<F> GateConstraint<F> for Logic<F>
where
    F: PrimeField,
{
    type CustomVals = LogicVals<F>;

    #[inline]
    fn constraints(
        separation_challenge: F,
        wit_vals: WitnessValues<F>,
        custom_vals: Self::CustomVals,
    ) -> F {
        let four = F::from(4_u64);
        let kappa = separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;
        let kappa_qu = kappa_cu * kappa;

        let a = custom_vals.a_next_eval - four * wit_vals.a_eval;
        let c_0 = delta(a);

        let b = custom_vals.b_next_eval - four * wit_vals.b_eval;
        let c_1 = delta(b) * kappa;

        let d = custom_vals.d_next_eval - four * wit_vals.d_eval;
        let c_2 = delta(d) * kappa_sq;

        let w = wit_vals.c_eval;
        let c_3 = (w - a * b) * kappa_cu;

        let c_4 = delta_xor_and(a, b, w, d, custom_vals.q_c_eval) * kappa_qu;

        (c_0 + c_1 + c_2 + c_3 + c_4) * separation_challenge
    }
}

/// Computes `f(f-1)(f-2)(f-3)`
pub(crate) fn delta<F>(f: F) -> F
where
    F: PrimeField,
{
    let f_1 = f - F::one();
    let f_2 = f - F::from(2_u64);
    let f_3 = f - F::from(3_u64);
    f * f_1 * f_2 * f_3
}

/// The identity we want to check is `q_logic * A = 0` where:
///
/// ```text
/// A = B + E
/// B = q_c * [9c - 3(a+b)]
/// E = 3(a+b+c) - 2F
/// F = w[w(4w - 18(a+b) + 81) + 18(a^2 + b^2) - 81(a+b) + 83]
/// ```
#[allow(non_snake_case)]
pub(crate) fn delta_xor_and<F>(a: F, b: F, w: F, c: F, q_c: F) -> F
where
    F: PrimeField,
{
    let nine = F::from(9_u64);
    let two = F::from(2_u64);
    let three = F::from(3_u64);
    let four = F::from(4_u64);
    let eighteen = F::from(18_u64);
    let eighty_one = F::from(81_u64);
    let eighty_three = F::from(83_u64);
    let F = w
        * (w * (four * w - eighteen * (a + b) + eighty_one)
            + eighteen * (a.square() + b.square())
            - eighty_one * (a + b)
            + eighty_three);
    let E = three * (a + b + c) - (two * F);
    let B = q_c * ((nine * c) - three * (a + b));
    B + E
}
