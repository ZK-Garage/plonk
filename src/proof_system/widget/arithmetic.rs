// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Arithmetic Gate

use crate::proof_system::linearisation_poly::ProofEvaluations;
use crate::proof_system::GateValues;
use ark_ec::PairingEngine;
use ark_ff::PrimeField;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::Evaluations;
use ark_poly::Polynomial;
use ark_poly_commit::sonic_pc::Commitment;
use ark_serialize::*;

/// Arithmetic Prover Key
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Eq, PartialEq)]
pub struct ProverKey<F>
where
    F: PrimeField,
{
    /// Arithmetic Selector
    pub arithmetic_selector: (DensePolynomial<F>, Evaluations<F>),

    /// Multiplication Selector
    pub mul_selector: (DensePolynomial<F>, Evaluations<F>),
}

impl<F> ProverKey<F>
where
    F: PrimeField,
{
    ///
    #[inline]
    fn constraint(mul_selector: F, values: GateValues<F>) -> F {
        (values.left * values.right * mul_selector)
            + (values.left * values.left_selector)
            + (values.right * values.right_selector)
            + (values.output * values.output_selector)
            + (values.fourth * values.fourth_selector)
            + values.constant_selector
    }

    ///
    #[inline]
    pub fn quotient_term(&self, index: usize, values: GateValues<F>) -> F {
        Self::constraint(self.mul_selector.1[index], values)
    }

    ///
    #[inline]
    pub fn linearisation_term(
        &self,
        z_challenge: &F,
        values: GateValues<F>,
    ) -> DensePolynomial<F> {
        &self.arithmetic_selector.0
            * Self::constraint(
                self.mul_selector.0.evaluate(z_challenge),
                values,
            )
    }
}

/// Arithmetic Verifier Key
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Eq, PartialEq)]
pub struct VerifierKey<E: PairingEngine> {
    /// Arithmetic Selector Commitment
    pub arithmetic_selector_commitment: Commitment<E>,

    /// Multiplication Selector Commitment
    pub mul_selector_commitment: Commitment<E>,
}

impl<E> VerifierKey<E>
where
    E: PairingEngine,
{
    ///
    pub fn extend_linearisation_commitment(
        &self,
        left_selector_commitment: &Commitment<E>,
        right_selector_commitment: &Commitment<E>,
        output_selector_commitment: &Commitment<E>,
        fourth_selector_commitment: &Commitment<E>,
        constant_selector_commitment: &Commitment<E>,
        evaluations: &ProofEvaluations<E::Fr>,
        scalars: &mut Vec<E::Fr>,
        points: &mut Vec<E::G1Affine>,
    ) {
        let q_arith_eval = evaluations.q_arith_eval;

        scalars.push(evaluations.a_eval * evaluations.b_eval * q_arith_eval);
        points.push(self.mul_selector_commitment.0);

        scalars.push(evaluations.a_eval * q_arith_eval);
        points.push(left_selector_commitment.0);

        scalars.push(evaluations.b_eval * q_arith_eval);
        points.push(right_selector_commitment.0);

        scalars.push(evaluations.c_eval * q_arith_eval);
        points.push(output_selector_commitment.0);

        scalars.push(evaluations.d_eval * q_arith_eval);
        points.push(fourth_selector_commitment.0);

        scalars.push(q_arith_eval);
        points.push(constant_selector_commitment.0);
    }
}
