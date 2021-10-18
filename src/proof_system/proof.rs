// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! A Proof stores the commitments to all of the elements that
//! are needed to univocally identify a prove of some statement.
//!
//! This module contains the implementation of the `StandardComposer`s
//! `Proof` structure and it's methods.

use super::{linearisation_poly::ProofEvaluations, PCAggregateProof};
use crate::commitment_scheme::kzg10::OpeningKey;
use crate::proof_system::VerifierKey;
use crate::{error::Error, transcript::TranscriptProtocol};
use ark_ec::{msm::VariableBaseMSM, AffineCurve, TEModelParameters};
use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_ff::{fields::batch_inversion, Field, FpParameters, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};
use ark_poly_commit::kzg10::Commitment;
use core::marker::PhantomData;
use merlin::Transcript;
use num_traits::{One, Zero};
use rayon::prelude::*;

/// A Proof is a composition of `Commitment`s to the Witness, Permutation,
/// Quotient, Shifted and Opening polynomials as well as the
/// `ProofEvaluations`.
///
/// It's main goal is to allow the `Verifier` to
/// formally verify that the secret witnesses used to generate the [`Proof`]
/// satisfy a circuit that both [`Prover`](super::Prover) and
/// [`Verifier`](super::Verifier) have in common succintly and without any
/// capabilities of adquiring any kind of knowledge about the witness used to
/// construct the Proof.
#[derive(Debug, Eq, PartialEq, Clone, Default)]
pub struct Proof<E: PairingEngine, P: TEModelParameters<BaseField = E::Fr>> {
    /// Commitment to the witness polynomial for the left wires.
    pub(crate) a_comm: Commitment<E>,
    /// Commitment to the witness polynomial for the right wires.
    pub(crate) b_comm: Commitment<E>,
    /// Commitment to the witness polynomial for the output wires.
    pub(crate) c_comm: Commitment<E>,
    /// Commitment to the witness polynomial for the fourth wires.
    pub(crate) d_comm: Commitment<E>,

    /// Commitment to the permutation polynomial.
    pub(crate) z_comm: Commitment<E>,

    /// Commitment to the quotient polynomial.
    pub(crate) t_1_comm: Commitment<E>,
    /// Commitment to the quotient polynomial.
    pub(crate) t_2_comm: Commitment<E>,
    /// Commitment to the quotient polynomial.
    pub(crate) t_3_comm: Commitment<E>,
    /// Commitment to the quotient polynomial.
    pub(crate) t_4_comm: Commitment<E>,

    /// Commitment to the opening polynomial.
    pub(crate) w_z_comm: Commitment<E>,
    /// Commitment to the shifted opening polynomial.
    pub(crate) w_zw_comm: Commitment<E>,
    /// Subset of all of the evaluations added to the proof.
    pub(crate) evaluations: ProofEvaluations<E::Fr>,
    _marker: PhantomData<P>,
}

impl<E: PairingEngine, P: TEModelParameters<BaseField = E::Fr>> Proof<E, P> {
    /// Performs the verification of a [`Proof`] returning a boolean result.
    pub(crate) fn verify(
        &self,
        verifier_key: &VerifierKey<E, P>,
        transcript: &mut Transcript,
        opening_key: &OpeningKey<E>,
        pub_inputs: &[E::Fr],
    ) -> Result<(), Error> {
        let domain =
            GeneralEvaluationDomain::<E::Fr>::new(verifier_key.n).unwrap();

        // Subgroup checks are done when the proof is deserialised.

        // In order for the Verifier and Prover to have the same view in the
        // non-interactive setting Both parties must commit the same
        // elements into the transcript Below the verifier will simulate
        // an interaction with the prover by adding the same elements
        // that the prover added into the transcript, hence generating the
        // same challenges
        //
        // Add commitment to witness polynomials to transcript
        transcript.append_commitment(b"w_l", &self.a_comm);
        transcript.append_commitment(b"w_r", &self.b_comm);
        transcript.append_commitment(b"w_o", &self.c_comm);
        transcript.append_commitment(b"w_4", &self.d_comm);

        // Compute beta and gamma challenges
        let beta = transcript.challenge_scalar(b"beta");
        transcript.append_scalar(b"beta", &beta);
        let gamma = transcript.challenge_scalar(b"gamma");
        // Add commitment to permutation polynomial to transcript
        transcript.append_commitment(b"z", &self.z_comm);

        // Compute quotient challenge
        let alpha = transcript.challenge_scalar(b"alpha");
        let range_sep_challenge =
            transcript.challenge_scalar(b"range separation challenge");
        let logic_sep_challenge =
            transcript.challenge_scalar(b"logic separation challenge");
        let fixed_base_sep_challenge =
            transcript.challenge_scalar(b"fixed base separation challenge");
        let var_base_sep_challenge =
            transcript.challenge_scalar(b"variable base separation challenge");

        // Add commitment to quotient polynomial to transcript
        transcript.append_commitment(b"t_1", &self.t_1_comm);
        transcript.append_commitment(b"t_2", &self.t_2_comm);
        transcript.append_commitment(b"t_3", &self.t_3_comm);
        transcript.append_commitment(b"t_4", &self.t_4_comm);

        // Compute evaluation challenge
        let z_challenge = transcript.challenge_scalar(b"z");

        // Compute zero polynomial evaluated at `z_challenge`
        let z_h_eval = domain.evaluate_vanishing_polynomial(z_challenge);

        // Compute first lagrange polynomial evaluated at `z_challenge`
        let l1_eval =
            compute_first_lagrange_evaluation(&domain, &z_h_eval, &z_challenge);

        // Compute quotient polynomial evaluated at `z_challenge`
        let t_eval = self.compute_quotient_evaluation(
            &domain,
            pub_inputs,
            alpha,
            beta,
            gamma,
            z_challenge,
            z_h_eval,
            l1_eval,
            self.evaluations.perm_eval,
        );

        // Compute commitment to quotient polynomial
        // This method is necessary as we pass the `un-splitted` variation
        // to our commitment scheme
        let t_comm =
            self.compute_quotient_commitment(&z_challenge, domain.size());

        // Add evaluations to transcript
        transcript.append_scalar(b"a_eval", &self.evaluations.a_eval);
        transcript.append_scalar(b"b_eval", &self.evaluations.b_eval);
        transcript.append_scalar(b"c_eval", &self.evaluations.c_eval);
        transcript.append_scalar(b"d_eval", &self.evaluations.d_eval);
        transcript.append_scalar(b"a_next_eval", &self.evaluations.a_next_eval);
        transcript.append_scalar(b"b_next_eval", &self.evaluations.b_next_eval);
        transcript.append_scalar(b"d_next_eval", &self.evaluations.d_next_eval);
        transcript
            .append_scalar(b"left_sig_eval", &self.evaluations.left_sigma_eval);
        transcript.append_scalar(
            b"right_sig_eval",
            &self.evaluations.right_sigma_eval,
        );
        transcript
            .append_scalar(b"out_sig_eval", &self.evaluations.out_sigma_eval);
        transcript
            .append_scalar(b"q_arith_eval", &self.evaluations.q_arith_eval);
        transcript.append_scalar(b"q_c_eval", &self.evaluations.q_c_eval);
        transcript.append_scalar(b"q_l_eval", &self.evaluations.q_l_eval);
        transcript.append_scalar(b"q_r_eval", &self.evaluations.q_r_eval);
        transcript.append_scalar(b"perm_eval", &self.evaluations.perm_eval);
        transcript.append_scalar(b"t_eval", &t_eval);
        transcript.append_scalar(b"r_eval", &self.evaluations.lin_poly_eval);

        // Compute linearisation commitment
        let r_comm = self.compute_linearisation_commitment(
            alpha,
            beta,
            gamma,
            (
                range_sep_challenge,
                logic_sep_challenge,
                fixed_base_sep_challenge,
                var_base_sep_challenge,
            ),
            z_challenge,
            l1_eval,
            &verifier_key,
        );

        // Commitment Scheme
        // Now we delegate computation to the commitment scheme by batch
        // checking two proofs.
        //
        // The `AggregateProof`, which proves that all the necessary
        // polynomials evaluated at `z_challenge` are
        // correct and a `Proof` which is proof that the
        // permutation polynomial evaluated at the shifted root of unity is
        // correct

        // Compose the Aggregated Proof
        //
        let mut aggregate_proof =
            PCAggregateProof::<E, P>::with_witness(self.w_z_comm);
        aggregate_proof.add_part((t_eval, t_comm));
        aggregate_proof.add_part((self.evaluations.lin_poly_eval, r_comm));
        aggregate_proof.add_part((self.evaluations.a_eval, self.a_comm));
        aggregate_proof.add_part((self.evaluations.b_eval, self.b_comm));
        aggregate_proof.add_part((self.evaluations.c_eval, self.c_comm));
        aggregate_proof.add_part((self.evaluations.d_eval, self.d_comm));
        aggregate_proof.add_part((
            self.evaluations.left_sigma_eval,
            verifier_key.permutation.left_sigma,
        ));
        aggregate_proof.add_part((
            self.evaluations.right_sigma_eval,
            verifier_key.permutation.right_sigma,
        ));
        aggregate_proof.add_part((
            self.evaluations.out_sigma_eval,
            verifier_key.permutation.out_sigma,
        ));
        // Flatten proof with opening challenge
        let flattened_proof_a = aggregate_proof.flatten(transcript);

        // Compose the shifted aggregate proof
        let mut shifted_aggregate_proof =
            PCAggregateProof::<E, P>::with_witness(self.w_zw_comm);
        shifted_aggregate_proof
            .add_part((self.evaluations.perm_eval, self.z_comm));
        shifted_aggregate_proof
            .add_part((self.evaluations.a_next_eval, self.a_comm));
        shifted_aggregate_proof
            .add_part((self.evaluations.b_next_eval, self.b_comm));
        shifted_aggregate_proof
            .add_part((self.evaluations.d_next_eval, self.d_comm));
        let flattened_proof_b = shifted_aggregate_proof.flatten(transcript);

        // Add commitment to openings to transcript
        transcript.append_commitment(b"w_z", &self.w_z_comm);
        transcript.append_commitment(b"w_z_w", &self.w_zw_comm);

        // XXX: Move this to a constants file with the K1...
        let group_gen = E::Fr::from_repr(
            <<E as PairingEngine>::Fr as PrimeField>::Params::GENERATOR,
        )
        .unwrap();

        // Batch check
        if opening_key
            .batch_check(
                &[z_challenge, (z_challenge * group_gen)],
                &[flattened_proof_a, flattened_proof_b],
                transcript,
            )
            .is_err()
        {
            return Err(Error::ProofVerificationError);
        }
        Ok(())
    }

    fn compute_quotient_evaluation(
        &self,
        domain: &GeneralEvaluationDomain<E::Fr>,
        pub_inputs: &[E::Fr],
        alpha: E::Fr,
        beta: E::Fr,
        gamma: E::Fr,
        z_challenge: E::Fr,
        z_h_eval: E::Fr,
        l1_eval: E::Fr,
        z_hat_eval: E::Fr,
    ) -> E::Fr {
        // Compute the public input polynomial evaluated at `z_challenge`
        let pi_eval = compute_barycentric_eval(pub_inputs, z_challenge, domain);

        let alpha_sq = alpha.square();
        // r + PI(z)
        let a = self.evaluations.lin_poly_eval + pi_eval;

        // a + beta * sigma_1 + gamma
        let beta_sig1 = beta * self.evaluations.left_sigma_eval;
        let b_0 = self.evaluations.a_eval + beta_sig1 + gamma;

        // b+ beta * sigma_2 + gamma
        let beta_sig2 = beta * self.evaluations.right_sigma_eval;
        let b_1 = self.evaluations.b_eval + beta_sig2 + gamma;

        // c+ beta * sigma_3 + gamma
        let beta_sig3 = beta * self.evaluations.out_sigma_eval;
        let b_2 = self.evaluations.c_eval + beta_sig3 + gamma;

        // ((d + gamma) * z_hat) * alpha
        let b_3 = (self.evaluations.d_eval + gamma) * z_hat_eval * alpha;

        let b = b_0 * b_1 * b_2 * b_3;

        // l_1(z) * alpha^2
        let c = l1_eval * alpha_sq;

        // Return t_eval
        (a - b - c) * z_h_eval.inverse().unwrap()
    }

    fn compute_quotient_commitment(
        &self,
        z_challenge: &E::Fr,
        n: usize,
    ) -> Commitment<E> {
        let z_n = z_challenge.pow(&[n as u64, 0, 0, 0]);
        let z_two_n = z_challenge.pow(&[2 * n as u64, 0, 0, 0]);
        let z_three_n = z_challenge.pow(&[3 * n as u64, 0, 0, 0]);
        let t_comm = self.t_1_comm.0.into_projective()
            + &(self.t_2_comm.0.mul(z_n.into_repr()))
            + &(self.t_3_comm.0.mul(z_two_n.into_repr()))
            + &(self.t_4_comm.0.mul(z_three_n.into_repr()));
        Commitment(t_comm.into_affine())
    }

    // Commitment to [r]_1
    fn compute_linearisation_commitment(
        &self,
        alpha: E::Fr,
        beta: E::Fr,
        gamma: E::Fr,
        (
            range_sep_challenge,
            logic_sep_challenge,
            fixed_base_sep_challenge,
            var_base_sep_challenge,
        ): (E::Fr, E::Fr, E::Fr, E::Fr),
        z_challenge: E::Fr,
        l1_eval: E::Fr,
        verifier_key: &VerifierKey<E, P>,
    ) -> Commitment<E> {
        let mut scalars: Vec<_> = Vec::with_capacity(6);
        let mut points: Vec<E::G1Affine> = Vec::with_capacity(6);

        verifier_key.arithmetic.compute_linearisation_commitment(
            &mut scalars,
            &mut points,
            &self.evaluations,
        );

        verifier_key.range.compute_linearisation_commitment(
            range_sep_challenge,
            &mut scalars,
            &mut points,
            &self.evaluations,
        );

        verifier_key.logic.compute_linearisation_commitment(
            logic_sep_challenge,
            &mut scalars,
            &mut points,
            &self.evaluations,
        );

        verifier_key.fixed_base.compute_linearisation_commitment(
            fixed_base_sep_challenge,
            &mut scalars,
            &mut points,
            &self.evaluations,
        );

        verifier_key.variable_base.compute_linearisation_commitment(
            var_base_sep_challenge,
            &mut scalars,
            &mut points,
            &self.evaluations,
        );

        verifier_key.permutation.compute_linearisation_commitment(
            &mut scalars,
            &mut points,
            &self.evaluations,
            z_challenge,
            (alpha, beta, gamma),
            l1_eval,
            self.z_comm.0,
        );

        let scalars_repr: Vec<<E::Fr as PrimeField>::BigInt> =
            scalars.iter().map(|s| s.into_repr()).collect();

        Commitment(
            VariableBaseMSM::multi_scalar_mul(&points, &scalars_repr).into(),
        )

        // Commitment::from_projective(VariableBaseMSM::multi_scalar_mul(
        //     &points,
        //     &scalars
        //         .into_iter()
        //         .map(|s| s.into_repr())
        //         .collect::<Vec<_>>(),
        // ))
    }
}

// TODO: Document this with the Lagrange formula if possible.
fn compute_first_lagrange_evaluation<F: PrimeField>(
    domain: &GeneralEvaluationDomain<F>,
    z_h_eval: &F,
    z_challenge: &F,
) -> F {
    let n_fr = F::from(domain.size() as u64);
    let denom = n_fr * &(*z_challenge - &F::one());
    *z_h_eval * denom.inverse().unwrap()
}

fn compute_barycentric_eval<F: PrimeField>(
    evaluations: &[F],
    point: F,
    domain: &GeneralEvaluationDomain<F>,
) -> F {
    let numerator = (point.pow(&[domain.size() as u64, 0, 0, 0]) - F::one())
        * F::from(domain.size() as u64).inverse().unwrap();
    let range = (0..evaluations.len()).into_par_iter();

    let non_zero_evaluations: Vec<usize> = range
        .filter(|&i| {
            let evaluation = &evaluations[i];
            evaluation != &F::zero()
        })
        .collect();

    // Only compute the denominators with non-zero evaluations
    let range = (0..non_zero_evaluations.len()).into_par_iter();

    let mut denominators: Vec<F> = range
        .clone()
        .map(|i| {
            // index of non-zero evaluation
            let index = non_zero_evaluations[i];

            // XXX: Move this to ctants file like K1...

            (F::from_repr(<F::Params>::GENERATOR)
                .unwrap()
                .inverse()
                .unwrap()
                .pow(&[index as u64, 0, 0, 0])
                * point)
                - F::one()
        })
        .collect();
    batch_inversion(&mut denominators);

    let result: F = range
        .map(|i| {
            let eval_index = non_zero_evaluations[i];
            let eval = evaluations[eval_index];

            denominators[i] * eval
        })
        .sum();

    result * numerator
}

/*
#[cfg(test)]
mod proof_tests {
    use super::*;
    use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
    use ark_poly_commit::kzg10::Commitment;
    use rand_core::OsRng;

    #[test]
    fn test_dusk_bytes_serde_proof() {
        let proof = Proof {
            a_comm: Commitment::<Bls12_381>::default(),
            b_comm: Commitment::<Bls12_381>::default(),
            c_comm: Commitment::<Bls12_381>::default(),
            d_comm: Commitment::<Bls12_381>::default(),
            z_comm: Commitment::<Bls12_381>::default(),
            t_1_comm: Commitment::<Bls12_381>::default(),
            t_2_comm: Commitment::<Bls12_381>::default(),
            t_3_comm: Commitment::<Bls12_381>::default(),
            t_4_comm: Commitment::<Bls12_381>::default(),
            w_z_comm: Commitment::<Bls12_381>::default(),
            w_zw_comm: Commitment::<Bls12_381>::default(),
            evaluations: ProofEvaluations {
                a_eval: BlsScalar::random(&mut OsRng),
                b_eval: BlsScalar::random(&mut OsRng),
                c_eval: BlsScalar::random(&mut OsRng),
                d_eval: BlsScalar::random(&mut OsRng),
                a_next_eval: BlsScalar::random(&mut OsRng),
                b_next_eval: BlsScalar::random(&mut OsRng),
                d_next_eval: BlsScalar::random(&mut OsRng),
                q_arith_eval: BlsScalar::random(&mut OsRng),
                q_c_eval: BlsScalar::random(&mut OsRng),
                q_l_eval: BlsScalar::random(&mut OsRng),
                q_r_eval: BlsScalar::random(&mut OsRng),
                left_sigma_eval: BlsScalar::random(&mut OsRng),
                right_sigma_eval: BlsScalar::random(&mut OsRng),
                out_sigma_eval: BlsScalar::random(&mut OsRng),
                lin_poly_eval: BlsScalar::random(&mut OsRng),
                perm_eval: BlsScalar::random(&mut OsRng),
            },
            _marker: PhantomData::new(),
        };

        let proof_bytes = proof.to_bytes();
        let got_proof = Proof::from_bytes(&proof_bytes).unwrap();
        assert_eq!(got_proof, proof);
    }
}
*/
