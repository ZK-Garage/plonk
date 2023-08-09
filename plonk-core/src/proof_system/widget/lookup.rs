// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) ZK-Garage. All rights reserved.
//! Lookup gates

use crate::error::Error;
use crate::lookup::multiset::MultiSet;
use crate::proof_system::linearisation_poly::ProofEvaluations;
use crate::util::lc;
use ark_ff::{FftField, PrimeField};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain};
use ark_poly_commit::PolynomialCommitment;
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
    pub table_1: MultiSet<F>,
    /// Column 2 of lookup table
    pub table_2: MultiSet<F>,
    /// Column 3 of lookup table
    pub table_3: MultiSet<F>,
    /// Column 4 of lookup table
    pub table_4: MultiSet<F>,
}

impl<F> ProverKey<F>
where
    F: PrimeField,
{
    /// Compute lookup portion of quotient polynomial
    pub fn compute_lookup_quotient_term(
        &self,
        domain: &GeneralEvaluationDomain<F>,
        wl_eval_8n: &[F],
        wr_eval_8n: &[F],
        wo_eval_8n: &[F],
        w4_eval_8n: &[F],
        f_eval_8n: &[F],
        table_eval_8n: &[F],
        h1_eval_8n: &[F],
        h2_eval_8n: &[F],
        z2_eval_8n: &[F],
        l1_eval_8n: &[F],
        delta: F,
        epsilon: F,
        zeta: F,
        lookup_sep: F,
    ) -> Result<Vec<F>, Error>
    where
        F: PrimeField,
    {
        let domain_8n = GeneralEvaluationDomain::<F>::new(8 * domain.size())
        .ok_or(Error::InvalidEvalDomainSize {
        log_size_of_group: (8 * domain.size()).trailing_zeros(),
        adicity:
            <<F as FftField>::FftParams as ark_ff::FftParameters>::TWO_ADICITY,
    })?;

        Ok((0..domain_8n.size())
            .map(|i| {
                self.compute_quotient_i(
                    i,
                    wl_eval_8n[i],
                    wr_eval_8n[i],
                    wo_eval_8n[i],
                    w4_eval_8n[i],
                    f_eval_8n[i],
                    table_eval_8n[i],
                    table_eval_8n[i + 8],
                    h1_eval_8n[i],
                    h1_eval_8n[i + 8],
                    h2_eval_8n[i],
                    z2_eval_8n[i],
                    z2_eval_8n[i + 8],
                    l1_eval_8n[i],
                    delta,
                    epsilon,
                    zeta,
                    lookup_sep,
                )
            })
            .collect())
    }

    /// Compute evals of lookup portion of quotient polynomial
    pub fn compute_quotient_i(
        &self,
        index: usize,
        w_l_i: F,
        w_r_i: F,
        w_o_i: F,
        w_4_i: F,
        f_i: F,
        table_i: F,
        table_i_next: F,
        h1_i: F,
        h1_i_next: F,
        h2_i: F,
        z2_i: F,
        z2_i_next: F,
        l1_i: F,
        delta: F,
        epsilon: F,
        zeta: F,
        lookup_sep: F,
    ) -> F {
        // q_lookup(X) * (a(X) + zeta * b(X) + (zeta^2 * c(X)) + (zeta^3 * d(X)
        // - f(X))) * α_1
        let lookup_sep_sq = lookup_sep.square();
        let lookup_sep_cu = lookup_sep_sq * lookup_sep;
        let one_plus_delta = delta + F::one();
        let epsilon_one_plus_delta = epsilon * one_plus_delta;

        let a = {
            let q_lookup_i = self.q_lookup.1[index];
            let compressed_tuple = lc(&[w_l_i, w_r_i, w_o_i, w_4_i], &zeta);
            q_lookup_i * (compressed_tuple - f_i) * lookup_sep
        };

        // z2(X) * (1+δ) * (ε+f(X)) * (ε*(1+δ) + t(X) + δt(Xω)) * lookup_sep^2
        let b = {
            let b_0 = epsilon + f_i;
            let b_1 = epsilon_one_plus_delta + table_i + delta * table_i_next;

            z2_i * one_plus_delta * b_0 * b_1 * lookup_sep_sq
        };

        // − z2(Xω) * (ε*(1+δ) + h1(X) + δ*h2(X)) * (ε*(1+δ) + h2(X) + δ*h1(Xω))
        // * lookup_sep^2
        let c = {
            let c_0 = epsilon_one_plus_delta + h1_i + delta * h2_i;
            let c_1 = epsilon_one_plus_delta + h2_i + delta * h1_i_next;

            -z2_i_next * c_0 * c_1 * lookup_sep_sq
        };

        let d = { (z2_i - F::one()) * l1_i * lookup_sep_cu };

        a + b + c + d
    }

    /// Compute linearization for lookup gates
    pub(crate) fn compute_linearisation(
        &self,
        l1_eval: F,
        a_eval: F,
        b_eval: F,
        c_eval: F,
        d_eval: F,
        f_eval: F,
        table_eval: F,
        table_next_eval: F,
        h1_next_eval: F,
        h2_eval: F,
        z2_next_eval: F,
        delta: F,
        epsilon: F,
        zeta: F,
        z2_poly: &DensePolynomial<F>,
        h1_poly: &DensePolynomial<F>,
        lookup_sep: F,
    ) -> DensePolynomial<F> {
        let lookup_sep_sq = lookup_sep.square();
        let lookup_sep_cu = lookup_sep * lookup_sep_sq;
        let one_plus_delta = delta + F::one();
        let epsilon_one_plus_delta = epsilon * one_plus_delta;

        let a = {
            let compressed_tuple = lc(&[a_eval, b_eval, c_eval, d_eval], &zeta);
            &self.q_lookup.0 * ((compressed_tuple - f_eval) * lookup_sep)
        };

        // z2(X) * (1 + δ) * (ε + f_bar) * (ε(1+δ) + t_bar + δ*tω_bar) *
        // lookup_sep^2
        let b = {
            let b_0 = epsilon + f_eval;
            let b_1 =
                epsilon_one_plus_delta + table_eval + delta * table_next_eval;
            let b_2 = l1_eval * lookup_sep_cu;

            z2_poly * (one_plus_delta * b_0 * b_1 * lookup_sep_sq + b_2)
        };

        // h1(X) * (−z2ω_bar) * (ε(1+δ) + h2_bar  + δh1ω_bar) * lookup_sep^2
        let c = {
            let c_0 = -z2_next_eval * lookup_sep_sq;
            let c_1 = epsilon_one_plus_delta + h2_eval + delta * h1_next_eval;

            h1_poly * (c_0 * c_1)
        };
        a + b + c
    }
}

/// LookUp Verifier Key
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone,
    Copy(bound = "PC::Commitment: Copy"),
    Debug(bound = "PC::Commitment: core::fmt::Debug"),
    Eq(bound = "PC::Commitment: Eq"),
    PartialEq(bound = "PC::Commitment: PartialEq")
)]
pub struct VerifierKey<F, PC>
where
    F: PrimeField,
    PC: PolynomialCommitment<F, DensePolynomial<F>>,
{
    /// Lookup Selector Commitment
    pub q_lookup: PC::Commitment,
    /// Commitment to first table column
    pub table_1: PC::Commitment,
    /// Commitment to second table column
    pub table_2: PC::Commitment,
    /// Commitment to third table column
    pub table_3: PC::Commitment,
    /// Commitment to fourth table column
    pub table_4: PC::Commitment,
}

impl<F, PC> VerifierKey<F, PC>
where
    F: PrimeField,
    PC: PolynomialCommitment<F, DensePolynomial<F>>,
{
    /// Computes the linearisation commitments.
    pub fn compute_linearisation_commitment(
        &self,
        scalars: &mut Vec<F>,
        points: &mut Vec<PC::Commitment>,
        evaluations: &ProofEvaluations<F>,
        (delta, epsilon, zeta): (F, F, F),
        lookup_sep: F,
        l1_eval: F,
        z2_comm: PC::Commitment,
        h1_comm: PC::Commitment,
    ) {
        let one_plus_delta = F::one() + delta;
        let epsilon_one_plus_delta = epsilon * one_plus_delta;
        let lookup_sep_sq = lookup_sep.square();
        let lookup_sep_cu = lookup_sep_sq * lookup_sep;

        let a = {
            let compressed_eval = lc(
                &[
                    evaluations.wire_evals.a_eval,
                    evaluations.wire_evals.b_eval,
                    evaluations.wire_evals.c_eval,
                    evaluations.wire_evals.d_eval,
                ],
                &zeta,
            );

            let a_0 = compressed_eval - evaluations.lookup_evals.f_eval;
            a_0 * lookup_sep
        };

        scalars.push(a);
        points.push(self.q_lookup.clone());

        // (1 + δ) * (ε + f_bar) * (ε(1+δ) + t_bar + δ*tω_bar) *  lookup_sep^2
        let b = {
            let b_0 = epsilon + evaluations.lookup_evals.f_eval;
            let b_1 = epsilon_one_plus_delta
                + evaluations.lookup_evals.table_eval
                + delta * evaluations.lookup_evals.table_next_eval;
            let b_2 = l1_eval * lookup_sep_cu;
            one_plus_delta * b_0 * b_1 * lookup_sep_sq + b_2
        };

        scalars.push(b);
        points.push(z2_comm);

        let c = {
            let c_0 = -evaluations.lookup_evals.z2_next_eval * lookup_sep_sq;
            let c_1 = epsilon_one_plus_delta
                + evaluations.lookup_evals.h2_eval
                + delta * evaluations.lookup_evals.h1_next_eval;
            c_0 * c_1
        };
        scalars.push(c);
        points.push(h1_comm);
    }
}
