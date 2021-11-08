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
    pub(crate) q_c: (DensePolynomial<F>, Evaluations<F>),
    pub(crate) q_logic: (DensePolynomial<F>, Evaluations<F>),
}

impl<F: PrimeField> ProverKey<F> {
    pub(crate) fn compute_quotient_i(
        &self,
        index: usize,
        logic_separation_challenge: F,
        w_l_i: F,
        w_l_i_next: F,
        w_r_i: F,
        w_r_i_next: F,
        w_o_i: F,
        w_4_i: F,
        w_4_i_next: F,
    ) -> F {
        let four = F::from(4 as u64);

        let q_logic_i = self.q_logic.1[index];
        let q_c_i = self.q_c.1[index];

        let kappa = logic_separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;
        let kappa_qu = kappa_cu * kappa;

        let a = w_l_i_next - four * w_l_i;
        let c_0 = delta(a);

        let b = w_r_i_next - four * w_r_i;
        let c_1 = delta(b) * kappa;

        let d = w_4_i_next - four * w_4_i;
        let c_2 = delta(d) * kappa_sq;

        let w = w_o_i;
        let c_3 = (w - a * b) * kappa_cu;

        let c_4 = delta_xor_and(a, b, w, d, q_c_i) * kappa_qu;

        q_logic_i * (c_3 + c_0 + c_1 + c_2 + c_4) * logic_separation_challenge
    }

    pub(crate) fn compute_linearisation(
        &self,
        logic_separation_challenge: F,
        a_eval: F,
        a_next_eval: F,
        b_eval: F,
        b_next_eval: F,
        c_eval: F,
        d_eval: F,
        d_next_eval: F,
        q_c_eval: F,
    ) -> DensePolynomial<F> {
        let four = F::from(4 as u64);
        let q_logic_poly = &self.q_logic.0;

        let kappa = logic_separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;
        let kappa_qu = kappa_cu * kappa;

        let a = a_next_eval - four * a_eval;
        let c_0 = delta(a);

        let b = b_next_eval - four * b_eval;
        let c_1 = delta(b) * kappa;

        let d = d_next_eval - four * d_eval;
        let c_2 = delta(d) * kappa_sq;

        let w = c_eval;
        let c_3 = (w - a * b) * kappa_cu;

        let c_4 = delta_xor_and(a, b, w, d, q_c_eval) * kappa_qu;

        let t = (c_0 + c_1 + c_2 + c_3 + c_4) * logic_separation_challenge;

        q_logic_poly * t
    }
}

// Computes f(f-1)(f-2)(f-3)
pub(crate) fn delta<F: PrimeField>(f: F) -> F {
    let f_1 = f - F::one();
    let f_2 = f - F::from(2 as u64);
    let f_3 = f - F::from(3 as u64);
    f * f_1 * f_2 * f_3
}

// The identity we want to check is q_logic * A = 0
// A = B + E
// B = q_c * [9c - 3(a+b)]
// E = 3(a+b+c) - 2F
// F = w[w(4w - 18(a+b) + 81) + 18(a^2 + b^2) - 81(a+b) + 83]
#[allow(non_snake_case)]
pub(crate) fn delta_xor_and<F: PrimeField>(
    a: F,
    b: F,
    w: F,
    c: F,
    q_c: F,
) -> F {
    let nine = F::from(9 as u64);
    let two = F::from(2 as u64);
    let three = F::from(3 as u64);
    let four = F::from(4 as u64);
    let eighteen = F::from(18 as u64);
    let eighty_one = F::from(81 as u64);
    let eighty_three = F::from(83 as u64);

    let F = w
        * (w * (four * w - eighteen * (a + b) + eighty_one)
            + eighteen * (a.square() + b.square())
            - eighty_one * (a + b)
            + eighty_three);
    let E = three * (a + b + c) - (two * F);
    let B = q_c * ((nine * c) - three * (a + b));
    B + E
}
