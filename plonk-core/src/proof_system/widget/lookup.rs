// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) ZK-Garage. All rights reserved.
//! Lookup gates

use crate::lookup::multiset::MultiSet;
use crate::proof_system::linearisation_poly::ProofEvaluations;
use crate::proof_system::widget::GateConstraint;
use crate::util::lc;
use ark_ff::PrimeField;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::Evaluations;
use ark_serialize::*;

/// Lookup Gates Prover Key
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Eq, PartialEq)]
pub struct ProverKey<F>
where
    F: PrimeField,
{
    /// Lookup selector
    pub q_lookup: (DensePolynomial<F>, Evaluations<F>),
    /// Column 1 of lookup table
    pub table_1: MultiSet<F>, // PC::Commitment, DensePolynomial<F>),
    /// Column 2 of lookup table
    pub table_2: MultiSet<F>, // PC::Commitment, DensePolynomial<F>),
    /// Column 3 of lookup table
    pub table_3: MultiSet<F>, // PC::Commitment, DensePolynomial<F>),
    /// Column 4 of lookup table
    pub table_4: MultiSet<F>, // PC::Commitment, DensePolynomial<F>),
}

impl<F> ProverKey<F>
where
    F: PrimeField,
{
    /// Compute lookup portion of quotient polynomial
    pub fn compute_quotient_i(
        &self,
        index: usize,
        w_l_i: F,
        w_r_i: F,
        w_o_i: F,
        w_4_i: F,
        f_i: F,
        zeta: F,
        lookup_challenge: F,
    ) -> F {
        // q_lookup(X) * (a(X) + zeta * b(X) + (zeta^2 * c(X)) + (zeta^3 * d(X)
        // - f(X))) * α_1
        let q_lookup_i = self.q_lookup.1[index];
        let compressed_tuple = Self::compress(w_l_i, w_r_i, w_o_i, w_4_i, zeta);

        q_lookup_i * (compressed_tuple - f_i) * lookup_challenge
    }

    /// Compute linearization for lookup gates
    pub(crate) fn compute_linearization(
        &self,
        a_eval: F,
        b_eval: F,
        c_eval: F,
        d_eval: F,
        f_eval: F,
        zeta: F,
        lookup_separation_challenge: F,
    ) -> DensePolynomial<F> {
        let l_sep_2 = lookup_separation_challenge.square();
        let l_sep_3 = l_sep_2 * lookup_separation_challenge;

        let a_0 = Self::compress(a_eval, b_eval, c_eval, d_eval, zeta);

        &self.q_lookup.0 * ((a_0 - f_eval) * l_sep_3)
    }

    /// Compresseses a row of values into a single field element
    /// by applying a random linear combination
    fn compress(w_l: F, w_r: F, w_o: F, w_4: F, zeta: F) -> F {
        lc(vec![w_l, w_r, w_o, w_4], zeta)
    }
}
