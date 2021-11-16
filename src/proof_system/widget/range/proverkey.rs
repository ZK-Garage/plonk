// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use ark_ff::PrimeField;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::Evaluations;
use ark_serialize::*;
#[derive(
    Debug, PartialEq, Eq, Clone, CanonicalDeserialize, CanonicalSerialize,
)]
pub(crate) struct ProverKey<F: PrimeField> {
    pub(crate) q_range: (DensePolynomial<F>, Evaluations<F>),
}

impl<F: PrimeField> ProverKey<F> {
    pub(crate) fn compute_quotient_i(
        &self,
        index: usize,
        range_separation_challenge: F,
        w_l_i: F,
        w_r_i: F,
        w_o_i: F,
        w_4_i: F,
        w_4_i_next: F,
    ) -> F {
        let four = F::from(4u64);
        let q_range_i = self.q_range.1[index];

        let kappa = range_separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;

        // Delta([c(X) - 4 * d(X)]) + Delta([b(X) - 4 * c(X)]) + Delta([a(X) - 4
        // * b(X)]) + Delta([d(Xg) - 4 * a(X)]) * Q_Range(X)
        //
        let b_1 = delta(w_o_i - four * w_4_i);
        let b_2 = delta(w_r_i - four * w_o_i) * kappa;
        let b_3 = delta(w_l_i - four * w_r_i) * kappa_sq;
        let b_4 = delta(w_4_i_next - four * w_l_i) * kappa_cu;
        (b_1 + b_2 + b_3 + b_4) * q_range_i * range_separation_challenge
    }

    pub(crate) fn compute_linearisation(
        &self,
        range_separation_challenge: F,
        a_eval: F,
        b_eval: F,
        c_eval: F,
        d_eval: F,
        d_next_eval: F,
    ) -> DensePolynomial<F> {
        let four = F::from(4u32);
        let q_range_poly = &self.q_range.0;

        let kappa = range_separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;

        // Delta([c_eval - 4 * d_eval]) + Delta([b_eval - 4 * c_eval]) +
        // Delta([a_eval - 4 * b_eval]) + Delta([d_next_eval - 4 * a_eval]) *
        // Q_Range(X)
        let b_1 = delta(c_eval - four * d_eval);
        let b_2 = delta(b_eval - four * c_eval) * kappa;
        let b_3 = delta(a_eval - four * b_eval) * kappa_sq;
        let b_4 = delta(d_next_eval - four * a_eval) * kappa_cu;

        let t = (b_1 + b_2 + b_3 + b_4) * range_separation_challenge;

        q_range_poly * t
    }
}

// Computes f(f-1)(f-2)(f-3)
pub(crate) fn delta<F: PrimeField>(f: F) -> F {
    let f_1 = f - F::one();
    let f_2 = f - F::from(2 as u64);
    let f_3 = f - F::from(3 as u64);
    f * f_1 * f_2 * f_3
}
