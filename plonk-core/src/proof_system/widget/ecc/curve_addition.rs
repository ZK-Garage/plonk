// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Elliptic Curve Point Addition Gate

use crate::proof_system::{
    widget::{GateConstraint, WitnessValues},
    CustomGateValues, ProverKey,
};
use ark_ec::{ModelParameters, TEModelParameters};
use ark_ff::{FftField, Field};
use ark_poly::{univariate::DensePolynomial, Polynomial};
use core::marker::PhantomData;
use std::collections::HashMap;

pub struct CAValues<F>
where
    F: Field,
{
    pub a_next: F,
    pub b_next: F,
    pub d_next: F,
}

impl<F> CustomGateValues<F> for CAValues<F>
where
    F: Field,
{
    fn new(vals: HashMap<String, F>) -> Self {
        let a_next = *vals.get(&"a_next".to_string()).unwrap();
        let b_next = *vals.get(&"b_next".to_string()).unwrap();
        let d_next = *vals.get(&"d_next".to_string()).unwrap();
        CAValues {
            a_next,
            b_next,
            d_next,
        }
    }
}

/// Curve Addition Gate
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct CurveAddition<F, P>(PhantomData<(F, P)>)
where
    F: FftField,
    P: ModelParameters<BaseField = F>;

impl<F, P> GateConstraint<F> for CurveAddition<F, P>
where
    F: FftField,
    P: TEModelParameters<BaseField = F>,
{
    type CustomValues = CAValues<F>;
    #[inline]
    fn constraints(
        separation_challenge: F,
        witness_vals: WitnessValues<F>,
        custom_vals: Self::CustomValues,
    ) -> F {
        let x_1 = witness_vals.left;
        let x_3 = custom_vals.a_next;
        let y_1 = witness_vals.right;
        let y_3 = custom_vals.b_next;
        let x_2 = witness_vals.output;
        let y_2 = witness_vals.fourth;
        let x1_y2 = custom_vals.d_next;

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

    fn evaluations(
        prover_key: &ProverKey<F>,
        w_l_poly: &DensePolynomial<F>,
        w_r_poly: &DensePolynomial<F>,
        w_o_poly: &DensePolynomial<F>,
        w_4_poly: &DensePolynomial<F>,
        z_challenge: &F,
        omega: F,
        custom_evals: HashMap<String, F>,
    ) {
        let shifted_z = *z_challenge * omega;

        let a_next_label = "a_next".to_string();
        let b_next_label = "b_next".to_string();
        let d_next_label = "d_next".to_string();

        if !custom_evals.contains_key(&a_next_label) {
            let a_next = w_l_poly.evaluate(&shifted_z);
            custom_evals.insert(a_next_label, a_next);
        }

        if !custom_evals.contains_key(&b_next_label) {
            let b_next = w_r_poly.evaluate(&shifted_z);
            custom_evals.insert(b_next_label, b_next);
        }

        if !custom_evals.contains_key(&d_next_label) {
            let d_next = w_4_poly.evaluate(&shifted_z);
            custom_evals.insert(d_next_label, d_next);
        }
    }
    fn verifier_key_term<PC>(
    ) -> Vec<ark_poly_commit::LabeledCommitment<PC::Commitment>>
    where
        PC: ark_poly_commit::PolynomialCommitment<
            F,
            ark_poly::univariate::DensePolynomial<F>,
        >,
    {
        unimplemented!();
    }
}
