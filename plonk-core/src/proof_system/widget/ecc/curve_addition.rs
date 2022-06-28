// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Elliptic Curve Point Addition Gate

use crate::{
    parameters::{CircuitParameters, EmbeddedCurve},
    proof_system::{
        ecc::{SWEmbeddedCurve, TEEmbeddedCurve},
        linearisation_poly::CustomEvaluations,
        widget::{GateConstraint, WitnessValues},
        CustomValues,
    },
};
use ark_ec::{SWModelParameters, TEModelParameters};
use ark_ff::PrimeField;
use core::marker::PhantomData;

/// Values needed for the computation of the Curve Addition gate constraint.
pub struct CAVals<F>
where
    F: PrimeField,
{
    /// Left wire value in the next position
    pub a_next_val: F,
    /// Right wire value in the next position
    pub b_next_val: F,
    /// Fourth wire value in the next position
    pub d_next_val: F,
}

impl<F> CustomValues<F> for CAVals<F>
where
    F: PrimeField,
{
    fn from_evaluations(custom_evals: &CustomEvaluations<F>) -> Self {
        let a_next_val = custom_evals.get("a_next_eval");
        let b_next_val = custom_evals.get("b_next_eval");
        let d_next_val = custom_evals.get("d_next_eval");
        CAVals {
            a_next_val,
            b_next_val,
            d_next_val,
        }
    }
}

/// Curve Addition Gate
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct CurveAddition<P, EC>(PhantomData<(P, EC)>)
where
    P: CircuitParameters,
    EC: EmbeddedCurve<P>;

impl<P, EmbeddedBaseField, EmbeddedCurveParameters>
    GateConstraint<P::ScalarField> for CurveAddition<P, TEEmbeddedCurve<P>>
where
    EmbeddedBaseField: PrimeField,
    EmbeddedCurveParameters: TEModelParameters<BaseField = EmbeddedBaseField>,
    P: CircuitParameters<
        ScalarField = EmbeddedBaseField,
        EmbeddedCurve = TEEmbeddedCurve<P>,
        EmbeddedCurveParameters = EmbeddedCurveParameters,
    >,
{
    type CustomVals = CAVals<P::ScalarField>;
    #[inline]
    fn constraints(
        separation_challenge: P::ScalarField,
        wit_vals: WitnessValues<P::ScalarField>,
        custom_vals: Self::CustomVals,
    ) -> P::ScalarField {
        let x_1 = wit_vals.a_val;
        let x_3 = custom_vals.a_next_val;
        let y_1 = wit_vals.b_val;
        let y_3 = custom_vals.b_next_val;
        let x_2 = wit_vals.c_val;
        let y_2 = wit_vals.d_val;
        let x1_y2 = custom_vals.d_next_val;

        let kappa = separation_challenge.square();

        // Check that `x1 * y2` is correct
        let xy_consistency = x_1 * y_2 - x1_y2;

        let y1_x2 = y_1 * x_2;
        let y1_y2 = y_1 * y_2;
        let x1_x2 = x_1 * x_2;

        // Check that `x_3` is correct
        let x3_lhs = x1_y2 + y1_x2;
        let x3_rhs =
            x_3 + (x_3 * EmbeddedCurveParameters::COEFF_D * x1_y2 * y1_x2);
        let x3_consistency = (x3_lhs - x3_rhs) * kappa;

        // Check that `y_3` is correct
        let y3_lhs = y1_y2 + x1_x2;
        let y3_rhs =
            y_3 - y_3 * EmbeddedCurveParameters::COEFF_D * x1_y2 * y1_x2;
        let y3_consistency = (y3_lhs - y3_rhs) * kappa.square();

        (xy_consistency + x3_consistency + y3_consistency)
            * separation_challenge
    }
}

impl<P, EmbeddedBaseField, EmbeddedCurveParameters>
    GateConstraint<P::ScalarField> for CurveAddition<P, SWEmbeddedCurve<P>>
where
    EmbeddedBaseField: PrimeField,
    EmbeddedCurveParameters: SWModelParameters<BaseField = EmbeddedBaseField>,
    P: CircuitParameters<
        ScalarField = EmbeddedBaseField,
        EmbeddedCurve = SWEmbeddedCurve<P>,
        EmbeddedCurveParameters = EmbeddedCurveParameters,
    >,
{
    type CustomVals = CAVals<P::ScalarField>;
    #[inline]
    fn constraints(
        separation_challenge: P::ScalarField,
        wit_vals: WitnessValues<P::ScalarField>,
        custom_vals: Self::CustomVals,
    ) -> P::ScalarField {
        let x_1 = wit_vals.a_val;
        let x_3 = custom_vals.a_next_val;
        let y_1 = wit_vals.b_val;
        let y_3 = custom_vals.b_next_val;
        let x_2 = wit_vals.c_val;
        let y_2 = wit_vals.d_val;
        //let x1_y2 = custom_vals.d_next_val;

        let kappa = separation_challenge.square();

        let x1_x2 = x_1 - x_2;
        let x2_x3 = x_2 - x_3;
        let y1_y2 = y_1 - y_2;
        let y2_y3 = y_2 + y_3;

        let constraint_1 = (x_1 + x_2 + x_3) * x1_x2 * x1_x2 - y1_y2 * y1_y2;
        let constraint_2 = x1_x2 * y2_y3 - x2_x3 * y1_y2;

        (constraint_1 + constraint_2 * kappa) * separation_challenge
    }
}
