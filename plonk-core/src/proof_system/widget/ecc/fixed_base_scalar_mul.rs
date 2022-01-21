// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Elliptic Curve Fixed-Base Scalar Multiplication Gate
//!
//! NOTE: The ECC gadget does not check that the initial point is on the
//! curve for two reasons:
//! - We constrain the accumulator to start from the identity point, which the
//!   verifier knows is on the curve
//! - We are adding multiples of the generator to the accumulator which the
//!   verifier also knows is on the curve and is prime order
//! - We do allow arbitrary BlsScalar multiplication, and possibly XXX: may add
//!   constraints to ensure the generator is correct (prime order)
//!
//! Bits are accumulated in base2. So we use d(Xw) - 2d(X) to extract the
//! base2 bit.

use crate::proof_system::{
    widget::{GateConstraint, WitnessValues},
    CustomGateValues, ProverKey,
};
use ark_ec::{ModelParameters, TEModelParameters};
use ark_ff::{FftField, Field};
use ark_poly::{univariate::DensePolynomial, Polynomial};
use core::marker::PhantomData;
use std::collections::HashMap;

pub struct FBSMValues<F>
where
    F: Field,
{
    pub a_next: F,
    pub b_next: F,
    pub d_next: F,
    pub q_r: F,
    pub q_l: F,
    pub q_c: F,
}

impl<F> CustomGateValues<F> for FBSMValues<F>
where
    F: Field,
{
    fn new(vals: HashMap<String, F>) -> Self {
        let a_next = *vals.get(&"a_next".to_string()).unwrap();
        let b_next = *vals.get(&"b_next".to_string()).unwrap();
        let d_next = *vals.get(&"d_next".to_string()).unwrap();
        let q_l = *vals.get(&"q_l".to_string()).unwrap();
        let q_r = *vals.get(&"q_r".to_string()).unwrap();
        let q_c = *vals.get(&"q_c".to_string()).unwrap();
        FBSMValues {
            a_next,
            b_next,
            d_next,
            q_l,
            q_r,
            q_c,
        }
    }
}

/// Fixed-Base Scalar Multiplication Gate
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct FixedBaseScalarMul<F, P>(PhantomData<(F, P)>)
where
    F: FftField,
    P: ModelParameters<BaseField = F>;

impl<F, P> GateConstraint<F> for FixedBaseScalarMul<F, P>
where
    F: FftField,
    P: TEModelParameters<BaseField = F>,
{
    type CustomValues = FBSMValues<F>;

    #[inline]
    fn constraints(
        separation_challenge: F,
        witness_vals: WitnessValues<F>,
        custom_vals: Self::CustomValues,
    ) -> F {
        let left_next = custom_vals.a_next;
        let right_next = custom_vals.b_next;
        let fourth_next = custom_vals.d_next;
        let left_selector = custom_vals.q_l;
        let right_selector = custom_vals.q_r;
        let constant_selector = custom_vals.q_r;

        let kappa = separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;

        let x_beta_eval = left_selector;
        let y_beta_eval = right_selector;

        let acc_x = witness_vals.left;
        let acc_x_next = left_next;
        let acc_y = witness_vals.right;
        let acc_y_next = right_next;

        let xy_alpha = witness_vals.output;

        let accumulated_bit = witness_vals.fourth;
        let accumulated_bit_next = fourth_next;
        let bit = extract_bit(accumulated_bit, accumulated_bit_next);

        // Check bit consistency
        let bit_consistency = check_bit_consistency(bit);

        let y_alpha = bit.square() * (y_beta_eval - F::one()) + F::one();
        let x_alpha = x_beta_eval * bit;

        // xy_alpha consistency check
        let xy_consistency = ((bit * constant_selector) - xy_alpha) * kappa;

        // x accumulator consistency check
        let x_3 = acc_x_next;
        let lhs = x_3 + (x_3 * xy_alpha * acc_x * acc_y * P::COEFF_D);
        let rhs = (x_alpha * acc_y) + (y_alpha * acc_x);
        let x_acc_consistency = (lhs - rhs) * kappa_sq;

        // y accumulator consistency check
        let y_3 = acc_y_next;
        let lhs = y_3 - (y_3 * xy_alpha * acc_x * acc_y * P::COEFF_D);
        let rhs = (x_alpha * acc_x) + (y_alpha * acc_y);
        let y_acc_consistency = (lhs - rhs) * kappa_cu;

        let checks = bit_consistency
            + x_acc_consistency
            + y_acc_consistency
            + xy_consistency;

        checks * separation_challenge
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
        let q_l_label = "q_l".to_string();
        let q_r_label = "q_r".to_string();
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

        if !custom_evals.contains_key(&q_l_label) {
            let q_l = prover_key.arithmetic.q_l.0.evaluate(z_challenge);
            custom_evals.insert(q_l_label, q_l);
        }

        if !custom_evals.contains_key(&q_r_label) {
            let q_r = prover_key.arithmetic.q_r.0.evaluate(z_challenge);
            custom_evals.insert(q_r_label, q_r);
        }

        if !custom_evals.contains_key(&q_c_label) {
            let q_c = prover_key.arithmetic.q_c.0.evaluate(z_challenge);
            custom_evals.insert(q_c_label, q_c);
        }
    }
}

/// Extracts the bit value from the accumulated bit.
pub(crate) fn extract_bit<F>(curr_acc: F, next_acc: F) -> F
where
    F: Field,
{
    next_acc - curr_acc - curr_acc
}

/// Ensures that the bit is either `+1`, `-1`, or `0`.
pub(crate) fn check_bit_consistency<F>(bit: F) -> F
where
    F: Field,
{
    let one = F::one();
    bit * (bit - one) * (bit + one)
}
