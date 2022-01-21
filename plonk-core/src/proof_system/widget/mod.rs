// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Proof System Widgets

pub mod arithmetic;
pub mod ecc;
pub mod logic;
pub mod range;

use crate::{
    commitment::HomomorphicCommitment,
    proof_system::{linearisation_poly::ProofEvaluations, permutation},
    transcript::TranscriptProtocol,
};
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::{univariate::DensePolynomial, Evaluations};
use ark_poly_commit::{
    LabeledCommitment, LabeledPolynomial, PolynomialCommitment,
};
use ark_serialize::*;
use std::collections::HashMap;

type LabeledDPolynomial<F> = LabeledPolynomial<F, DensePolynomial<F>>;

/// Witness Values
///
/// This data structures holds the wire values for a given gate.
#[derive(derivative::Derivative)]
#[derivative(Clone, Debug, Default, Eq, Hash, PartialEq)]
pub struct WitnessValues<F>
where
    F: Field,
{
    /// Left Value
    pub left: F,

    /// Right Value
    pub right: F,

    /// Output Value
    pub output: F,

    /// Fourth Value
    pub fourth: F,
}

/// Witness Values
///
/// Data structures that hold the extra necessary values for a custom gate
pub trait CustomGateValues<F>
where
    F: Field,
{
    fn new(vals: HashMap<String, F>) -> Self;
}

/// Gate Constraint
pub trait GateConstraint<F>
where
    F: FftField,
{
    type CustomValues: CustomGateValues<F>;
    /// Returns the coefficient of the quotient polynomial for this gate given
    /// an instantiation of the gate at `values` and a
    /// `separation_challenge` if this gate requires it for soundness.
    ///
    /// This method is an encoding of the polynomial constraint(s) that this
    /// gate represents whenever it is added to a circuit.
    fn constraints(
        separation_challenge: F,
        witness_vals: WitnessValues<F>,
        custom_vals: Self::CustomValues,
    ) -> F;

    /// Computes the quotient polynomial term for the given gate type for the
    /// given value of `selector` instantiated with `separation_challenge` and
    /// `values`.
    fn quotient_term(
        selector: F,
        separation_challenge: F,
        witness_vals: WitnessValues<F>,
        custom_vals: Self::CustomValues,
    ) -> F {
        selector
            * Self::constraints(separation_challenge, witness_vals, custom_vals)
    }

    /// Computes the linearisation polynomial term for the given gate type
    /// at the `selector_polynomial` instantiated with `separation_challenge`
    /// and `values`.
    fn linearisation_term(
        selector_polynomial: &DensePolynomial<F>,
        separation_challenge: F,
        witness_vals: WitnessValues<F>,
        custom_vals: Self::CustomValues,
    ) -> DensePolynomial<F> {
        selector_polynomial
            * Self::constraints(separation_challenge, witness_vals, custom_vals)
    }

    /// Computes the necessary evaluations associated to the custom gate.
    /// These evaluations are the apportation of the custom gate to the proof.
    fn evaluations(
        prover_key: &ProverKey<F>,
        w_l_poly: &DensePolynomial<F>,
        w_r_poly: &DensePolynomial<F>,
        w_o_poly: &DensePolynomial<F>,
        w_4_poly: &DensePolynomial<F>,
        z_challenge: &F,
        omega: F,
        custom_evals: HashMap<String, F>,
    );

    /// Computes the commitments for the verifier key associated to the given
    /// gate type
    fn verifier_key_term<PC>() -> Vec<LabeledCommitment<PC::Commitment>>
    where
        PC: PolynomialCommitment<F, DensePolynomial<F>>,
    {
        unimplemented!();
    }

    /// Extends `scalars` and `points` to build the linearisation commitment
    /// with the given instantiation of `evaluations` and
    /// `separation_challenge`.
    fn extend_linearisation_commitment<PC>(
        selector_commitment: &PC::Commitment,
        separation_challenge: F,
        evaluations: &ProofEvaluations<F>,
        scalars: &mut Vec<F>,
        points: &mut Vec<PC::Commitment>,
    ) where
        PC: PolynomialCommitment<F, DensePolynomial<F>>,
    {
        let coefficient = Self::constraints(
            separation_challenge,
            WitnessValues {
                left: evaluations.a_eval,
                right: evaluations.b_eval,
                output: evaluations.c_eval,
                fourth: evaluations.d_eval,
            },
        );
        scalars.push(coefficient);
        points.push(selector_commitment.clone());
    }
}

/// PLONK circuit Verification Key.
///
/// This structure is used by the Verifier in order to verify a
/// [`Proof`](super::Proof).
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(
        bound = "arithmetic::VerifierKey<F,PC>: std::fmt::Debug, PC::Commitment: std::fmt::Debug"
    ),
    Eq(bound = "arithmetic::VerifierKey<F,PC>: Eq, PC::Commitment: Eq"),
    PartialEq(
        bound = "arithmetic::VerifierKey<F,PC>: PartialEq, PC::Commitment: PartialEq"
    )
)]
pub struct VerifierKey<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    /// Circuit size (not padded to a power of two).
    pub(crate) n: usize,

    /// Arithmetic Verifier Key
    pub(crate) arithmetic: arithmetic::VerifierKey<F, PC>,

    /// VerifierKey for permutation checks
    pub(crate) permutation: permutation::VerifierKey<PC::Commitment>,

    /// Labeled commitment to custom gates polynomials
    pub(crate) custom_gate: Vec<LabeledCommitment<PC::Commitment>>,
}

impl<F, PC> VerifierKey<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    /// Constructs a [`VerifierKey`] from the widget VerifierKey's that are
    /// constructed based on the selector polynomial commitments and the
    /// sigma polynomial commitments.
    pub(crate) fn from_polynomial_commitments(
        n: usize,
        q_m: PC::Commitment,
        q_l: PC::Commitment,
        q_r: PC::Commitment,
        q_o: PC::Commitment,
        q_4: PC::Commitment,
        q_c: PC::Commitment,
        q_arith: PC::Commitment,

        left_sigma: PC::Commitment,
        right_sigma: PC::Commitment,
        out_sigma: PC::Commitment,
        fourth_sigma: PC::Commitment,

        custom_gates: Vec<LabeledCommitment<PC::Commitment>>,
    ) -> Self {
        Self {
            n,
            arithmetic: arithmetic::VerifierKey {
                q_m,
                q_l,
                q_r,
                q_o,
                q_4,
                q_c,
                q_arith,
            },
            permutation: permutation::VerifierKey {
                left_sigma,
                right_sigma,
                out_sigma,
                fourth_sigma,
            },
            custom_gate: custom_gates,
        }
    }

    /// Returns the Circuit size padded to the next power of two.
    pub fn padded_circuit_size(&self) -> usize {
        self.n.next_power_of_two()
    }
}

impl<F, PC> VerifierKey<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    /// Adds the circuit description to the transcript.
    pub(crate) fn seed_transcript<T>(&self, transcript: &mut T)
    where
        T: TranscriptProtocol,
    {
        transcript.append(b"q_m", &self.arithmetic.q_m);
        transcript.append(b"q_l", &self.arithmetic.q_l);
        transcript.append(b"q_r", &self.arithmetic.q_r);
        transcript.append(b"q_o", &self.arithmetic.q_o);
        transcript.append(b"q_c", &self.arithmetic.q_c);
        transcript.append(b"q_4", &self.arithmetic.q_4);

        transcript.append(b"q_arith", &self.arithmetic.q_arith);
        transcript.append(b"left_sigma", &self.permutation.left_sigma);
        transcript.append(b"right_sigma", &self.permutation.right_sigma);
        transcript.append(b"out_sigma", &self.permutation.out_sigma);
        transcript.append(b"fourth_sigma", &self.permutation.fourth_sigma);

        self.custom_gate.iter().map(|labeled_commitment| {
            transcript.append(
                labeled_commitment.label().clone().as_bytes(),
                labeled_commitment.commitment(),
            );
        });

        transcript.circuit_domain_sep(self.n as u64);
    }
}

/// PLONK circuit Proving Key.
///
/// This structure is used by the Prover in order to construct a
/// [`Proof`](crate::proof_system::Proof).
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
    Eq(bound = ""),
    PartialEq(bound = "")
)]
pub struct ProverKey<F>
where
    F: FftField,
{
    /// Circuit size
    pub(crate) n: usize,

    /// Arithmetic Prover Key
    pub(crate) arithmetic: arithmetic::ProverKey<F>,

    /// ProverKey for permutation checks
    pub(crate) permutation: permutation::ProverKey<F>,

    pub(crate) custom_gates:
        HashMap<String, (DensePolynomial<F>, Evaluations<F>)>,

    /// Pre-processes the 8n Evaluations for the vanishing polynomial, so
    /// they do not need to be computed at the proving stage.
    ///
    /// NOTE: With this, we can combine all parts of the quotient polynomial
    /// in their evaluation phase and divide by the quotient
    /// polynomial without having to perform IFFT
    pub(crate) v_h_coset_8n: Evaluations<F>,
}

impl<F> ProverKey<F>
where
    F: FftField,
{
    pub(crate) fn v_h_coset_8n(&self) -> &Evaluations<F> {
        &self.v_h_coset_8n
    }

    /// Constructs a [`ProverKey`] from the widget ProverKey's that are
    /// constructed based on the selector polynomials and the
    /// sigma polynomials and it's evaluations.
    pub(crate) fn from_polynomials_and_evals(
        n: usize,
        q_m: (DensePolynomial<F>, Evaluations<F>),
        q_l: (DensePolynomial<F>, Evaluations<F>),
        q_r: (DensePolynomial<F>, Evaluations<F>),
        q_o: (DensePolynomial<F>, Evaluations<F>),
        q_4: (DensePolynomial<F>, Evaluations<F>),
        q_c: (DensePolynomial<F>, Evaluations<F>),
        q_arith: (DensePolynomial<F>, Evaluations<F>),
        left_sigma: (DensePolynomial<F>, Evaluations<F>),
        right_sigma: (DensePolynomial<F>, Evaluations<F>),
        out_sigma: (DensePolynomial<F>, Evaluations<F>),
        fourth_sigma: (DensePolynomial<F>, Evaluations<F>),

        custom_gates: HashMap<String, (DensePolynomial<F>, Evaluations<F>)>,

        linear_evaluations: Evaluations<F>,
        v_h_coset_8n: Evaluations<F>,
    ) -> Self {
        Self {
            n,
            arithmetic: arithmetic::ProverKey {
                q_m,
                q_l,
                q_r,
                q_o,
                q_4,
                q_c,
                q_arith,
            },
            permutation: permutation::ProverKey {
                left_sigma,
                right_sigma,
                out_sigma,
                fourth_sigma,
                linear_evaluations,
            },
            custom_gates,
            v_h_coset_8n,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::batch_test;
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_ec::models::TEModelParameters;
    use ark_poly::polynomial::univariate::DensePolynomial;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, UVPolynomial};
    use rand::rngs::OsRng;

    fn rand_labeled_poly_eval<F>(
        n: usize,
        label: &str,
    ) -> (LabeledDPolynomial<F>, Evaluations<F>)
    where
        F: FftField,
    {
        let polynomial = LabeledPolynomial::new(
            label.to_string(),
            DensePolynomial::rand(n, &mut OsRng),
            None,
            None,
        );
        (polynomial, rand_evaluations(n))
    }

    fn rand_poly_eval<F>(n: usize) -> (DensePolynomial<F>, Evaluations<F>)
    where
        F: FftField,
    {
        let polynomial = DensePolynomial::rand(n, &mut OsRng);
        (polynomial, rand_evaluations(n))
    }

    fn rand_evaluations<F>(n: usize) -> Evaluations<F>
    where
        F: FftField,
    {
        let domain = GeneralEvaluationDomain::new(8 * n).unwrap();
        let values: Vec<_> = (0..8 * n).map(|_| F::rand(&mut OsRng)).collect();
        Evaluations::from_vec_and_domain(values, domain)
    }

    #[test]
    fn test_serialise_deserialise_prover_key() {
        type F = ark_bls12_381::Fr;
        let n = 1 << 11;

        let q_m = rand_poly_eval(n);
        let q_l = rand_poly_eval(n);
        let q_r = rand_poly_eval(n);
        let q_o = rand_poly_eval(n);
        let q_4 = rand_poly_eval(n);
        let q_c = rand_poly_eval(n);
        let q_arith = rand_poly_eval(n);

        let left_sigma = rand_poly_eval(n);
        let right_sigma = rand_poly_eval(n);
        let out_sigma = rand_poly_eval(n);
        let fourth_sigma = rand_poly_eval(n);

        let custom_gates = vec![
            rand_labeled_poly_eval(n, "q_range"),
            rand_labeled_poly_eval(n, "q_logic"),
            rand_labeled_poly_eval(n, "q_fixed_group_add"),
            rand_labeled_poly_eval(n, "q_variable_group_add"),
        ];

        let linear_evaluations = rand_evaluations(n);
        let v_h_coset_8n = rand_evaluations(n);

        let prover_key = ProverKey::from_polynomials_and_evals(
            n,
            q_m,
            q_l,
            q_r,
            q_o,
            q_4,
            q_c,
            q_arith,
            left_sigma,
            right_sigma,
            out_sigma,
            fourth_sigma,
            custom_gates,
            linear_evaluations,
            v_h_coset_8n,
        );

        let mut prover_key_bytes = vec![];
        prover_key
            .serialize_unchecked(&mut prover_key_bytes)
            .unwrap();

        let obtained_pk: ProverKey<F> =
            ProverKey::deserialize_unchecked(prover_key_bytes.as_slice())
                .unwrap();

        assert_eq!(prover_key, obtained_pk);
    }

    fn test_serialise_deserialise_verifier_key<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
        VerifierKey<F, PC>: PartialEq,
    {
        let n = 2usize.pow(5);

        let q_m = PC::Commitment::default();
        let q_l = PC::Commitment::default();
        let q_r = PC::Commitment::default();
        let q_o = PC::Commitment::default();
        let q_4 = PC::Commitment::default();
        let q_c = PC::Commitment::default();
        let q_arith = PC::Commitment::default();

        let left_sigma = PC::Commitment::default();
        let right_sigma = PC::Commitment::default();
        let out_sigma = PC::Commitment::default();
        let fourth_sigma = PC::Commitment::default();

        let custom_gates = vec![
            LabeledCommitment::new(
                "q_range".to_string(),
                PC::Commitment::default(),
                None,
            ),
            LabeledCommitment::new(
                "q_logic".to_string(),
                PC::Commitment::default(),
                None,
            ),
            LabeledCommitment::new(
                "q_fixed_group_add".to_string(),
                PC::Commitment::default(),
                None,
            ),
            LabeledCommitment::new(
                "q_variable_group_add".to_string(),
                PC::Commitment::default(),
                None,
            ),
        ];

        let verifier_key = VerifierKey::<F, PC>::from_polynomial_commitments(
            n,
            q_m,
            q_l,
            q_r,
            q_o,
            q_4,
            q_c,
            q_arith,
            left_sigma,
            right_sigma,
            out_sigma,
            fourth_sigma,
            custom_gates,
        );

        let mut verifier_key_bytes = vec![];
        verifier_key
            .serialize_unchecked(&mut verifier_key_bytes)
            .unwrap();

        let obtained_vk: VerifierKey<F, PC> =
            VerifierKey::deserialize_unchecked(verifier_key_bytes.as_slice())
                .unwrap();

        assert!(verifier_key == obtained_vk);
    }

    // Test for Bls12_381
    batch_test!(
        [test_serialise_deserialise_verifier_key],
        [] => (
            Bls12_381, ark_ed_on_bls12_381::EdwardsParameters      )
    );

    // Test for Bls12_377
    batch_test!(
        [test_serialise_deserialise_verifier_key],
        [] => (
            Bls12_377, ark_ed_on_bls12_377::EdwardsParameters       )
    );
}
