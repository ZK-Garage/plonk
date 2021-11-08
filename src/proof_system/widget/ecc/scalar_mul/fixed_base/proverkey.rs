// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use ark_ec::TEModelParameters;
use ark_ff::PrimeField;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::Evaluations;
use ark_serialize::*;
use core::marker::PhantomData;
#[derive(
    Debug, PartialEq, Eq, Clone, CanonicalDeserialize, CanonicalSerialize,
)]
pub(crate) struct ProverKey<F: PrimeField, P: TEModelParameters<BaseField = F>>
{
    pub(crate) q_l: (DensePolynomial<F>, Evaluations<F>),
    pub(crate) q_r: (DensePolynomial<F>, Evaluations<F>),
    pub(crate) q_c: (DensePolynomial<F>, Evaluations<F>),
    pub(crate) q_fixed_group_add: (DensePolynomial<F>, Evaluations<F>),
    pub(crate) _marker: PhantomData<P>,
}

impl<F: PrimeField, P: TEModelParameters<BaseField = F>> ProverKey<F, P> {
    pub(crate) fn compute_quotient_i(
        &self,
        index: usize,
        ecc_separation_challenge: F,
        w_l_i: F,      // acc_x or curr_x
        w_l_i_next: F, //  // next_x
        w_r_i: F,      // acc_y or curr_y
        w_r_i_next: F, // next_y
        w_o_i: F,      // xy_alpha
        w_4_i: F,      // accumulated_bit
        w_4_i_next: F, // accumulated_bit_next
    ) -> F {
        let q_fixed_group_add_i = &self.q_fixed_group_add.1[index];
        let q_c_i = &self.q_c.1[index];

        let kappa = ecc_separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;

        let x_beta = &self.q_l.1[index];
        let y_beta = &self.q_r.1[index];

        let acc_x = w_l_i;
        let acc_x_next = w_l_i_next;
        let acc_y = w_r_i;
        let acc_y_next = w_r_i_next;

        let xy_alpha = w_o_i;

        let accumulated_bit = w_4_i;
        let accumulated_bit_next = w_4_i_next;
        let bit = extract_bit(accumulated_bit, accumulated_bit_next);

        // Checks
        //
        // Check bit consistency
        let bit_consistency = check_bit_consistency(bit);

        // Derive y_alpha and x_alpha from bit
        let y_alpha = bit.square() * (*y_beta - F::one()) + F::one();
        let x_alpha = bit * x_beta;

        // xy_alpha consistency check
        let xy_consistency = ((bit * q_c_i) - xy_alpha) * kappa;

        // x accumulator consistency check
        let x_3 = acc_x_next;
        let lhs = x_3 + (x_3 * xy_alpha * acc_x * acc_y * P::COEFF_D);
        let rhs = (acc_x * y_alpha) + (acc_y * x_alpha);
        let x_acc_consistency = (lhs - rhs) * kappa_sq;

        // y accumulator consistency check
        let y_3 = acc_y_next;
        let lhs = y_3 - (y_3 * xy_alpha * acc_x * acc_y * P::COEFF_D);
        let rhs = (acc_y * y_alpha) + (acc_x * x_alpha);
        let y_acc_consistency = (lhs - rhs) * kappa_cu;

        let identity = bit_consistency
            + x_acc_consistency
            + y_acc_consistency
            + xy_consistency;

        identity * q_fixed_group_add_i * ecc_separation_challenge
    }

    pub(crate) fn compute_linearisation(
        &self,
        ecc_separation_challenge: F,
        a_eval: F,
        a_next_eval: F,
        b_eval: F,
        b_next_eval: F,
        c_eval: F,
        d_eval: F,
        d_next_eval: F,
        q_l_eval: F,
        q_r_eval: F,
        q_c_eval: F,
    ) -> DensePolynomial<F> {
        let q_fixed_group_add_poly = &self.q_fixed_group_add.0;

        let kappa = ecc_separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;

        let x_beta_eval = q_l_eval;
        let y_beta_eval = q_r_eval;

        let acc_x = a_eval;
        let acc_x_next = a_next_eval;
        let acc_y = b_eval;
        let acc_y_next = b_next_eval;

        let xy_alpha = c_eval;

        let accumulated_bit = d_eval;
        let accumulated_bit_next = d_next_eval;
        let bit = extract_bit(accumulated_bit, accumulated_bit_next);

        // Check bit consistency
        let bit_consistency = check_bit_consistency(bit);

        let y_alpha = bit.square() * (y_beta_eval - F::one()) + F::one();

        let x_alpha = x_beta_eval * bit;

        // xy_alpha consistency check
        let xy_consistency = ((bit * q_c_eval) - xy_alpha) * kappa;

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

        let a = bit_consistency
            + x_acc_consistency
            + y_acc_consistency
            + xy_consistency;

        q_fixed_group_add_poly * (a * ecc_separation_challenge)
    }
}

pub(crate) fn extract_bit<F: PrimeField>(curr_acc: F, next_acc: F) -> F {
    // Next - 2 * current
    next_acc - curr_acc - curr_acc
}

// Ensures that the bit is either +1, -1 or 0
pub(crate) fn check_bit_consistency<F: PrimeField>(bit: F) -> F {
    let one = F::one();
    bit * (bit - one) * (bit + one)
}
