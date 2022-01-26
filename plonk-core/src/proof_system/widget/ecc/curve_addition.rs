// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Elliptic Curve Point Addition Gate

use crate::{
    get_label,
    proof_system::{
        linearisation_poly::CustomEvaluations,
        widget::{GateConstraint, WitnessValues},
        CustomValues,
    },
};
use ark_ec::{ModelParameters, TEModelParameters};
use ark_ff::PrimeField;
use core::marker::PhantomData;

pub struct CAVals<F>
where
    F: PrimeField,
{
    pub a_next_eval: F,
    pub b_next_eval: F,
    pub d_next_eval: F,
}

impl<F> CustomValues<F> for CAVals<F>
where
    F: PrimeField,
{
    fn from_evaluations(custom_evals: CustomEvaluations<F>) -> Self {
        let a_next_eval = custom_evals.get(get_label!(a_next_eval));
        let b_next_eval = custom_evals.get(get_label!(b_next_eval));
        let d_next_eval = custom_evals.get(get_label!(d_next_eval));
        CAVals {
            a_next_eval,
            b_next_eval,
            d_next_eval,
        }
    }
}

/// Curve Addition Gate
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct CurveAddition<F, P>(PhantomData<(F, P)>)
where
    F: PrimeField,
    P: ModelParameters<BaseField = F>;

impl<F, P> GateConstraint<F> for CurveAddition<F, P>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    type CustomVals = CAVals<F>;
    #[inline]
    fn constraints(
        separation_challenge: F,
        wit_vals: WitnessValues<F>,
        custom_vals: Self::CustomVals,
    ) -> F {
        let x_1 = wit_vals.a_eval;
        let x_3 = custom_vals.a_next_eval;
        let y_1 = wit_vals.r_eval;
        let y_3 = custom_vals.b_next_eval;
        let x_2 = wit_vals.c_eval;
        let y_2 = wit_vals.d_eval;
        let x1_y2 = custom_vals.d_next_eval;

        let kappa = separation_challenge.square();

        // Check that `x1 * y2` is correct
        let xy_consistency = x_1 * y_2 - x1_y2;

        let y1_x2 = y_1 * x_2;
        let y1_y2 = y_1 * y_2;
        let x1_x2 = x_1 * x_2;

        // Check that `x_3` is correct
        let x3_lhs = x1_y2 + y1_x2;
        let x3_rhs = x_3 + (x_3 * P::COEFF_D * x1_y2 * y1_x2);
        let x3_consistency = (x3_lhs - x3_rhs) * kappa;

        // Check that `y_3` is correct
        let y3_lhs = y1_y2 + x1_x2;
        let y3_rhs = y_3 - y_3 * P::COEFF_D * x1_y2 * y1_x2;
        let y3_consistency = (y3_lhs - y3_rhs) * kappa.square();

        (xy_consistency + x3_consistency + y3_consistency)
            * separation_challenge
    }
}
