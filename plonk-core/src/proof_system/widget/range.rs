// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Range Gate

use crate::proof_system::{GateConstraint, WitnessValues};
use ark_ff::{FftField, Field};
use ark_poly::{univariate::DensePolynomial, Polynomial};
use core::marker::PhantomData;
use std::collections::HashMap;

use super::{CustomGateValues, ProverKey};

pub struct RangeValues<F>
where
    F: Field,
{
    pub d_next: F,
}

impl<F> CustomGateValues<F> for RangeValues<F>
where
    F: Field,
{
    fn new(vals: HashMap<String, F>) -> Self {
        let d_next = *vals.get(&"d_next".to_string()).unwrap();
        RangeValues { d_next }
    }
}

/// Range Gate
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct Range<F>(PhantomData<F>)
where
    F: FftField;

impl<F> GateConstraint<F> for Range<F>
where
    F: FftField,
{
    type CustomValues = RangeValues<F>;
    #[inline]
    fn constraints(
        separation_challenge: F,
        witness_vals: WitnessValues<F>,
        custom_vals: Self::CustomValues,
    ) -> F {
        let four = F::from(4u64);
        let kappa = separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;

        //TODO Handle errors
        let b_1 = delta(witness_vals.output - four * witness_vals.fourth);
        let b_2 =
            delta(witness_vals.right - four * witness_vals.output) * kappa;
        let b_3 =
            delta(witness_vals.left - four * witness_vals.right) * kappa_sq;
        let b_4 =
            delta(custom_vals.d_next - four * witness_vals.left) * kappa_cu;
        (b_1 + b_2 + b_3 + b_4) * separation_challenge
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

        let d_next_label = "d_next".to_string();

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

/// Computes `f(f-1)(f-2)(f-3)`.
fn delta<F>(f: F) -> F
where
    F: Field,
{
    let f_1 = f - F::one();
    let f_2 = f - F::from(2_u64);
    let f_3 = f - F::from(3_u64);
    f * f_1 * f_2 * f_3
}
