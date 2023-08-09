// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Arithmetic Gates

use crate::proof_system::WitnessValues;
use crate::{
    constraint_system::SBOX_ALPHA,
    proof_system::linearisation_poly::ProofEvaluations,
};
use ark_ff::{FftField, PrimeField};
use ark_poly::{polynomial::univariate::DensePolynomial, Evaluations};
use ark_poly_commit::PolynomialCommitment;
use ark_serialize::*;

/// Arithmetic Gates Prover Key
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Eq, PartialEq)]
pub struct ProverKey<F>
where
    F: FftField,
{
    /// Multiplication Selector
    pub q_m: (DensePolynomial<F>, Evaluations<F>),

    /// Left Wire Selector
    pub q_l: (DensePolynomial<F>, Evaluations<F>),

    /// Right Wire Selector
    pub q_r: (DensePolynomial<F>, Evaluations<F>),

    /// Output Wire Selector
    pub q_o: (DensePolynomial<F>, Evaluations<F>),

    /// Fourth Wire Selector
    pub q_4: (DensePolynomial<F>, Evaluations<F>),

    /// Constant Selector
    pub q_c: (DensePolynomial<F>, Evaluations<F>),

    /// High degree selector
    pub q_hl: (DensePolynomial<F>, Evaluations<F>),

    /// High degree selector
    pub q_hr: (DensePolynomial<F>, Evaluations<F>),

    /// High degree selector
    pub q_h4: (DensePolynomial<F>, Evaluations<F>),

    /// Arithmetic Selector
    pub q_arith: (DensePolynomial<F>, Evaluations<F>),
}

impl<F> ProverKey<F>
where
    F: PrimeField,
{
    /// Computes the arithmetic gate contribution to the quotient polynomial at
    /// the element of the domain at the given `index`.
    pub fn compute_quotient_i(
        &self,
        index: usize,
        wit_vals: WitnessValues<F>,
    ) -> F {
        ((wit_vals.a_val * wit_vals.b_val * self.q_m.1[index])
            + (wit_vals.a_val * self.q_l.1[index])
            + (wit_vals.b_val * self.q_r.1[index])
            + (wit_vals.c_val * self.q_o.1[index])
            + (wit_vals.d_val * self.q_4.1[index])
            + (wit_vals.a_val.pow([SBOX_ALPHA]) * self.q_hl.1[index])
            + (wit_vals.b_val.pow([SBOX_ALPHA]) * self.q_hr.1[index])
            + (wit_vals.d_val.pow([SBOX_ALPHA]) * self.q_h4.1[index])
            + self.q_c.1[index])
            * self.q_arith.1[index]
    }

    /// Computes the arithmetic gate contribution to the linearisation
    /// polynomial at the given evaluation points.
    pub fn compute_linearisation(
        &self,
        a_eval: F,
        b_eval: F,
        c_eval: F,
        d_eval: F,
        q_arith_eval: F,
    ) -> DensePolynomial<F> {
        &(&((&self.q_m.0 * (a_eval * b_eval))
            + (&self.q_l.0 * a_eval)
            + (&self.q_r.0 * b_eval)
            + (&self.q_o.0 * c_eval)
            + (&self.q_4.0 * d_eval)
            + (&self.q_hl.0 * a_eval.pow([SBOX_ALPHA]))
            + (&self.q_hr.0 * b_eval.pow([SBOX_ALPHA]))
            + (&self.q_h4.0 * d_eval.pow([SBOX_ALPHA])))
            + &self.q_c.0)
            * q_arith_eval
    }
}

/// Arithmetic Gates Verifier Key
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone,
    Copy(bound = "PC::Commitment: Copy"),
    Debug(bound = "PC::Commitment: std::fmt::Debug"),
    Eq(bound = "PC::Commitment: Eq"),
    PartialEq(bound = "PC::Commitment: PartialEq")
)]
pub struct VerifierKey<F, PC>
where
    F: PrimeField,
    PC: PolynomialCommitment<F, DensePolynomial<F>>,
{
    /// Multiplication Selector Commitment
    pub q_m: PC::Commitment,

    /// Left Selector Commitment
    pub q_l: PC::Commitment,

    /// Right Selector Commitment
    pub q_r: PC::Commitment,

    /// Output Selector Commitment
    pub q_o: PC::Commitment,

    /// Fourth Selector Commitment
    pub q_4: PC::Commitment,

    /// Constant Selector Commitment
    pub q_c: PC::Commitment,

    /// High degree left Selector Commitment
    pub q_hl: PC::Commitment,

    /// High degree right Selector Commitment
    pub q_hr: PC::Commitment,

    /// High degree 4-th Selector Commitment
    pub q_h4: PC::Commitment,

    /// Arithmetic Selector Commitment
    pub q_arith: PC::Commitment,
}

impl<F, PC> VerifierKey<F, PC>
where
    F: PrimeField,
    PC: PolynomialCommitment<F, DensePolynomial<F>>,
{
    /// Computes arithmetic gate contribution to the linearisation polynomial
    /// commitment.
    pub fn compute_linearisation_commitment(
        &self,
        scalars: &mut Vec<F>,
        points: &mut Vec<PC::Commitment>,
        evaluations: &ProofEvaluations<F>,
    ) {
        let q_arith_eval = evaluations.custom_evals.get("q_arith_eval");

        scalars.push(
            evaluations.wire_evals.a_eval
                * evaluations.wire_evals.b_eval
                * q_arith_eval,
        );
        points.push(self.q_m.clone());

        scalars.push(evaluations.wire_evals.a_eval * q_arith_eval);
        points.push(self.q_l.clone());

        scalars.push(evaluations.wire_evals.b_eval * q_arith_eval);
        points.push(self.q_r.clone());

        scalars.push(evaluations.wire_evals.d_eval * q_arith_eval);
        points.push(self.q_4.clone());

        scalars.push(evaluations.wire_evals.c_eval * q_arith_eval);
        points.push(self.q_o.clone());

        scalars.push(
            evaluations.wire_evals.a_eval.pow([SBOX_ALPHA]) * q_arith_eval,
        );
        points.push(self.q_hl.clone());

        scalars.push(
            evaluations.wire_evals.b_eval.pow([SBOX_ALPHA]) * q_arith_eval,
        );
        points.push(self.q_hr.clone());

        scalars.push(
            evaluations.wire_evals.d_eval.pow([SBOX_ALPHA]) * q_arith_eval,
        );
        points.push(self.q_h4.clone());

        scalars.push(q_arith_eval);
        points.push(self.q_c.clone());
    }
}
