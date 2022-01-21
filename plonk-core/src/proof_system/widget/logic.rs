// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved..

//! Logic Gates

use crate::proof_system::widget::{GateConstraint, WitnessValues};
use ark_ff::{FftField, Field};
use ark_poly::{univariate::DensePolynomial, Polynomial};
use ark_poly_commit::{LabeledCommitment, PolynomialCommitment};
use core::marker::PhantomData;
use std::collections::HashMap;

use super::{CustomGateValues, ProverKey};

pub struct LogicValues<F>
where
    F: Field,
{
    pub a_next: F,
    pub b_next: F,
    pub d_next: F,
    pub q_c: F,
}

impl<F> CustomGateValues<F> for LogicValues<F>
where
    F: Field,
{
    fn new(vals: HashMap<String, F>) -> Self {
        let a_next = *vals.get(&"a_next".to_string()).unwrap();
        let b_next = *vals.get(&"b_next".to_string()).unwrap();
        let d_next = *vals.get(&"d_next".to_string()).unwrap();
        let q_c = *vals.get(&"q_c".to_string()).unwrap();
        LogicValues {
            a_next,
            b_next,
            d_next,
            q_c,
        }
    }
}

/// Logic Gate
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct Logic<F>(PhantomData<F>)
where
    F: FftField;

impl<F> GateConstraint<F> for Logic<F>
where
    F: FftField,
{
    type CustomValues = LogicValues<F>;
    #[inline]
    fn constraints(
        separation_challenge: F,
        witness_vals: WitnessValues<F>,
        custom_vals: Self::CustomValues,
    ) -> F {
        let four = F::from(4_u64);
        let kappa = separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;
        let kappa_qu = kappa_cu * kappa;

        let a = custom_vals.a_next - four * witness_vals.left;
        let c_0 = delta(a);

        let b = custom_vals.b_next - four * witness_vals.right;
        let c_1 = delta(b) * kappa;

        let d = custom_vals.d_next - four * witness_vals.fourth;
        let c_2 = delta(d) * kappa_sq;

        let w = witness_vals.output;
        let c_3 = (w - a * b) * kappa_cu;

        let c_4 = delta_xor_and(a, b, w, d, custom_vals.q_c) * kappa_qu;

        (c_0 + c_1 + c_2 + c_3 + c_4) * separation_challenge
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
        let q_c_label = "q_c".to_string();

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

        if !custom_evals.contains_key(&q_c_label) {
            let q_c = prover_key.arithmetic.q_c.0.evaluate(&z_challenge);
            custom_evals.insert(q_c_label, q_c);
        }
    }
    fn verifier_key_term<PC>() -> Vec<LabeledCommitment<PC::Commitment>>
    where
        PC: PolynomialCommitment<F, DensePolynomial<F>>,
    {
        unimplemented!();
    }
}

/// Computes `f(f-1)(f-2)(f-3)`
pub(crate) fn delta<F>(f: F) -> F
where
    F: Field,
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
    F: Field,
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
