// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) ZK-Garage. All rights reserved.

use crate::lookup::multiset::MultiSet;
use crate::proof_system::linearisation_poly::ProofEvaluations;
use crate::proof_system::widget::GateConstraint;
use crate::proof_system::widget::HomomorphicCommitment;
use crate::proof_system::WitnessValues;
use ark_ec::PairingEngine;
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::Evaluations;
use ark_serialize::*;
use core::marker::PhantomData;

/// Lookup Gates Prover Key
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = "PC::Commitment: std::fmt::Debug"),
    Eq(bound = "PC::Commitment: Eq"),
    PartialEq(bound = "PC::Commitment: PartialEq")
)]

/// Lookup gates prover key
pub struct ProverKey<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    pub q_lookup: (DensePolynomial<F>, Evaluations<F>),
    pub table_1: (MultiSet<F>, PC::Commitment, DensePolynomial<F>),
    pub table_2: (MultiSet<F>, PC::Commitment, DensePolynomial<F>),
    pub table_3: (MultiSet<F>, PC::Commitment, DensePolynomial<F>),
    pub table_4: (MultiSet<F>, PC::Commitment, DensePolynomial<F>),
}

impl<F, PC> ProverKey<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    pub fn compute_quotient_i(
        &self,
        index: usize,
        wit_vals: WitnessValues<F>,
        f_i: F,
        zeta: F,
        lookup_challenge: F,
    ) -> F {
        // q_lookup(X) * (a(X) + zeta * b(X) + (zeta^2 * c(X)) + (zeta^3 * d(X)
        // - f(X))) * α_1
        let a = {
            let q_lookup_i = self.q_lookup.1[index];
            let compressed_tuple = Self::compress(
                wit_vals.a_val,
                wit_vals.b_val,
                wit_vals.c_val,
                wit_vals.d_val,
                zeta,
            );

            q_lookup_i * (compressed_tuple - f_i) * lookup_challenge
        };

        a
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
        let zeta_sq = zeta * zeta;
        let zeta_cu = zeta * zeta_sq;

        // q_lookup(X) * f_eval * lookup_separation_challenge
        let a = {
            let a_0 =
                a_eval + zeta * b_eval + zeta_sq * c_eval + zeta_cu * d_eval;

            &self.q_lookup.0 * ((a_0 - f_eval) * l_sep_3)
        };

        a
    }

    fn compress(w_l: F, w_r: F, w_o: F, w_4: F, zeta: F) -> F {
        let zeta_sq = zeta.square();
        let zeta_cu = zeta_sq * zeta;

        let a = w_l;

        let b = w_r * zeta;

        let c = w_o * zeta_sq;

        let d = w_4 * zeta_cu;

        a + b + c + d
    }
}
