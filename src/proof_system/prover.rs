// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Prover-side of the PLONK Proving System

use crate::{
    constraint_system::{StandardComposer, Variable},
    error::Error,
    proof_system::{
        linearisation_poly, proof::Proof, quotient_poly, ProverKey,
    },
    transcript::TranscriptProtocol,
    util,
};
use ark_ec::{PairingEngine, TEModelParameters};
use ark_ff::{FftField, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
    UVPolynomial,
};
use ark_poly_commit::kzg10::{Powers, KZG10};
use core::ops::Add;
use merlin::Transcript;
use num_traits::Zero;

/// Abstraction structure designed to construct a circuit and generate
/// [`Proof`]s for it.
pub struct Prover<F>
where
    F: FftField,
{
    /// Proving Key which is used to create proofs about a specific PLONK
    /// circuit.
    pub prover_key: Option<ProverKey<F>>,

    /// Circuit Description
    pub(crate) cs: StandardComposer<F>,

    /// Store the messages exchanged during the preprocessing stage.
    ///
    /// This is copied each time, we make a proof.
    pub preprocessed_transcript: Transcript,
}

impl<F> Prover<F>
where
    F: FftField,
{
    /// Creates a new `Prover` instance.
    pub fn new(label: &'static [u8]) -> Self {
        Self {
            prover_key: None,
            cs: StandardComposer::new(),
            preprocessed_transcript: Transcript::new(label),
        }
    }

    /// Creates a new `Prover` object with some expected size.
    pub fn with_expected_size(label: &'static [u8], size: usize) -> Self {
        Self {
            prover_key: None,
            cs: StandardComposer::with_expected_size(size),
            preprocessed_transcript: Transcript::new(label),
        }
    }

    /// Returns a mutable copy of the underlying [`StandardComposer`].
    pub fn mut_cs(&mut self) -> &mut StandardComposer<F> {
        &mut self.cs
    }

    /// Returns the number of gates in the circuit thet the `Prover` actually
    /// stores inside.
    pub fn circuit_size(&self) -> usize {
        self.cs.circuit_size()
    }

    /// Preprocesses the underlying constraint system.
    pub fn preprocess<E: PairingEngine<Fr = F>>(
        &mut self,
        commit_key: &Powers<E>,
    ) -> Result<(), Error>
    where
        F: PrimeField,
    {
        if self.prover_key.is_some() {
            return Err(Error::CircuitAlreadyPreprocessed);
        }
        let pk = self
            .cs
            .preprocess_prover(commit_key, &mut self.preprocessed_transcript)?;
        self.prover_key = Some(pk);
        Ok(())
    }

    /// Split `t(X)` poly into 4 degree `n` polynomials.
    #[allow(clippy::type_complexity)] // NOTE: This is an ok type for internal use.
    fn split_tx_poly(
        &self,
        n: usize,
        t_x: &DensePolynomial<F>,
    ) -> (
        DensePolynomial<F>,
        DensePolynomial<F>,
        DensePolynomial<F>,
        DensePolynomial<F>,
    ) {
        (
            DensePolynomial::from_coefficients_vec(t_x[0..n].to_vec()),
            DensePolynomial::from_coefficients_vec(t_x[n..2 * n].to_vec()),
            DensePolynomial::from_coefficients_vec(t_x[2 * n..3 * n].to_vec()),
            DensePolynomial::from_coefficients_vec(t_x[3 * n..].to_vec()),
        )
    }

    /// Computes the quotient Opening [`DensePolynomial`].
    fn compute_quotient_opening_poly(
        n: usize,
        t_1_poly: &DensePolynomial<F>,
        t_2_poly: &DensePolynomial<F>,
        t_3_poly: &DensePolynomial<F>,
        t_4_poly: &DensePolynomial<F>,
        z_challenge: &F,
    ) -> DensePolynomial<F> {
        // Compute z^n , z^2n , z^3n
        let z_n = z_challenge.pow(&[n as u64, 0, 0, 0]);
        let z_two_n = z_challenge.pow(&[2 * n as u64, 0, 0, 0]);
        let z_three_n = z_challenge.pow(&[3 * n as u64, 0, 0, 0]);
        let a = t_1_poly;
        let b = t_2_poly * z_n;
        let c = t_3_poly * z_two_n;
        let d = t_4_poly * z_three_n;
        a + &b + c + d
    }

    /// Convert variables to their actual witness values.
    fn to_scalars(&self, vars: &[Variable]) -> Vec<F> {
        vars.iter().map(|var| self.cs.variables[var]).collect()
    }

    /// Resets the witnesses in the prover object.
    ///
    /// This function is used when the user wants to make multiple proofs with
    /// the same circuit.
    pub fn clear_witness(&mut self) {
        self.cs = StandardComposer::new();
    }

    /// Clears all data in the `Prover` instance.
    ///
    /// This function is used when the user wants to use the same `Prover` to
    /// make a [`Proof`] regarding a different circuit.
    pub fn clear(&mut self) {
        self.clear_witness();
        self.prover_key = None;
        self.preprocessed_transcript = Transcript::new(b"plonk");
    }

    /// Keys the [`Transcript`] with additional seed information
    /// Wrapper around [`Transcript::append_message`].
    ///
    /// [`Transcript`]: merlin::Transcript
    /// [`Transcript::append_message`]: merlin::Transcript::append_message
    pub fn key_transcript(&mut self, label: &'static [u8], message: &[u8]) {
        self.preprocessed_transcript.append_message(label, message);
    }

    /// Computes a single witness for multiple polynomials at the same point, by
    /// taking a random linear combination of the individual witnesses.
    ///
    /// The result does not depend on `z`, thus we can remove the term `f(z)`.
    fn compute_aggregate_witness(
        polynomials: &[DensePolynomial<F>],
        point: &F,
        challenge: F,
    ) -> DensePolynomial<F> {
        util::ruffini(
            util::powers_of(challenge)
                .zip(polynomials)
                .map(|(challenge, poly)| poly * challenge)
                .fold(Zero::zero(), Add::add),
            *point,
        )
    }

    /// Creates a [`Proof]` that demonstrates that a circuit is satisfied.
    /// # Note
    /// If you intend to construct multiple [`Proof`]s with different witnesses,
    /// after calling this method, the user should then call
    /// [`Prover::clear_witness`].
    /// This is automatically done when [`Prover::prove`] is called.
    pub fn prove_with_preprocessed<E: PairingEngine<Fr = F>, P>(
        &self,
        commit_key: &Powers<E>,
        prover_key: &ProverKey<F>,
    ) -> Result<Proof<E>, Error>
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
    {
        let domain =
            GeneralEvaluationDomain::new(self.cs.circuit_size()).unwrap();

        // Since the caller is passing a pre-processed circuit
        // We assume that the Transcript has been seeded with the preprocessed
        // Commitments
        let mut transcript = self.preprocessed_transcript.clone();

        // 1. Compute witness Polynomials
        //
        // Convert Variables to scalars padding them to the
        // correct domain size.
        let pad = vec![F::zero(); domain.size() - self.cs.w_l.len()];
        let w_l_scalar = &[&self.to_scalars(&self.cs.w_l)[..], &pad].concat();
        let w_r_scalar = &[&self.to_scalars(&self.cs.w_r)[..], &pad].concat();
        let w_o_scalar = &[&self.to_scalars(&self.cs.w_o)[..], &pad].concat();
        let w_4_scalar = &[&self.to_scalars(&self.cs.w_4)[..], &pad].concat();

        // Witnesses are now in evaluation form, convert them to coefficients
        // so that we may commit to them.
        let w_l_poly =
            DensePolynomial::from_coefficients_vec(domain.ifft(w_l_scalar));
        let w_r_poly =
            DensePolynomial::from_coefficients_vec(domain.ifft(w_r_scalar));
        let w_o_poly =
            DensePolynomial::from_coefficients_vec(domain.ifft(w_o_scalar));
        let w_4_poly =
            DensePolynomial::from_coefficients_vec(domain.ifft(w_4_scalar));

        // Commit to witness polynomials.
        let w_l_poly_commit = KZG10::commit(commit_key, &w_l_poly, None, None)?;
        let w_r_poly_commit = KZG10::commit(commit_key, &w_r_poly, None, None)?;
        let w_o_poly_commit = KZG10::commit(commit_key, &w_o_poly, None, None)?;
        let w_4_poly_commit = KZG10::commit(commit_key, &w_4_poly, None, None)?;

        // Add witness polynomial commitments to transcript.
        transcript.append(b"w_l", &w_l_poly_commit.0);
        transcript.append(b"w_r", &w_r_poly_commit.0);
        transcript.append(b"w_o", &w_o_poly_commit.0);
        transcript.append(b"w_4", &w_4_poly_commit.0);

        // 2. Compute permutation polynomial
        //
        // Compute permutation challenges; `beta` and `gamma`.
        let beta = transcript.challenge_scalar(b"beta");
        transcript.append(b"beta", &beta);
        let gamma = transcript.challenge_scalar(b"gamma");

        assert!(beta != gamma, "challenges must be different");

        let z_poly = DensePolynomial::from_coefficients_slice(
            &self.cs.perm.compute_permutation_poly(
                &domain,
                (w_l_scalar, w_r_scalar, w_o_scalar, w_4_scalar),
                beta,
                gamma,
                (
                    &prover_key.permutation.left_sigma.0,
                    &prover_key.permutation.right_sigma.0,
                    &prover_key.permutation.out_sigma.0,
                    &prover_key.permutation.fourth_sigma.0,
                ),
            ),
        );

        // Commit to permutation polynomial.
        let z_poly_commit = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key, &z_poly, None, None,
        )?;

        // Add permutation polynomial commitment to transcript.
        transcript.append(b"z", &z_poly_commit.0);

        // 3. Compute public inputs polynomial.
        let pi_poly = DensePolynomial::from_coefficients_vec(
            domain.ifft(&self.cs.construct_dense_pi_vec()),
        );

        // 4. Compute quotient polynomial
        //
        // Compute quotient challenge; `alpha`, and gate-specific separation
        // challenges.
        let alpha = transcript.challenge_scalar(b"alpha");
        let range_sep_challenge =
            transcript.challenge_scalar(b"range separation challenge");
        let logic_sep_challenge =
            transcript.challenge_scalar(b"logic separation challenge");
        let fixed_base_sep_challenge =
            transcript.challenge_scalar(b"fixed base separation challenge");
        let var_base_sep_challenge =
            transcript.challenge_scalar(b"variable base separation challenge");

        let t_poly = quotient_poly::compute::<F, P>(
            &domain,
            prover_key,
            &z_poly,
            &w_l_poly,
            &w_r_poly,
            &w_o_poly,
            &w_4_poly,
            &pi_poly,
            &alpha,
            &beta,
            &gamma,
            &range_sep_challenge,
            &logic_sep_challenge,
            &fixed_base_sep_challenge,
            &var_base_sep_challenge,
        )?;

        // Split quotient polynomial into 4 degree `n` polynomials
        let (t_1_poly, t_2_poly, t_3_poly, t_4_poly) =
            self.split_tx_poly(domain.size(), &t_poly);

        // Commit to splitted quotient polynomial
        let t_1_commit = KZG10::commit(commit_key, &t_1_poly, None, None)?;
        let t_2_commit = KZG10::commit(commit_key, &t_2_poly, None, None)?;
        let t_3_commit = KZG10::commit(commit_key, &t_3_poly, None, None)?;
        let t_4_commit = KZG10::commit(commit_key, &t_4_poly, None, None)?;

        // Add quotient polynomial commitments to transcript
        transcript.append(b"t_1", &t_1_commit.0);
        transcript.append(b"t_2", &t_2_commit.0);
        transcript.append(b"t_3", &t_3_commit.0);
        transcript.append(b"t_4", &t_4_commit.0);

        // 4. Compute linearisation polynomial
        //
        // Compute evaluation challenge; `z`.
        let z_challenge = transcript.challenge_scalar(b"z");

        let (lin_poly, evaluations) = linearisation_poly::compute::<F, P>(
            &domain,
            prover_key,
            &alpha,
            &beta,
            &gamma,
            &range_sep_challenge,
            &logic_sep_challenge,
            &fixed_base_sep_challenge,
            &var_base_sep_challenge,
            &z_challenge,
            &w_l_poly,
            &w_r_poly,
            &w_o_poly,
            &w_4_poly,
            &t_poly,
            &z_poly,
        );

        // Add evaluations to transcript.
        transcript.append(b"a_eval", &evaluations.proof.a_eval);
        transcript.append(b"b_eval", &evaluations.proof.b_eval);
        transcript.append(b"c_eval", &evaluations.proof.c_eval);
        transcript.append(b"d_eval", &evaluations.proof.d_eval);
        transcript.append(b"a_next_eval", &evaluations.proof.a_next_eval);
        transcript.append(b"b_next_eval", &evaluations.proof.b_next_eval);
        transcript.append(b"d_next_eval", &evaluations.proof.d_next_eval);
        transcript.append(b"left_sig_eval", &evaluations.proof.left_sigma_eval);
        transcript
            .append(b"right_sig_eval", &evaluations.proof.right_sigma_eval);
        transcript.append(b"out_sig_eval", &evaluations.proof.out_sigma_eval);
        transcript.append(b"q_arith_eval", &evaluations.proof.q_arith_eval);
        transcript.append(b"q_c_eval", &evaluations.proof.q_c_eval);
        transcript.append(b"q_l_eval", &evaluations.proof.q_l_eval);
        transcript.append(b"q_r_eval", &evaluations.proof.q_r_eval);
        transcript.append(b"perm_eval", &evaluations.proof.permutation_eval);
        transcript.append(b"t_eval", &evaluations.quot_eval);
        transcript.append(
            b"r_eval",
            &evaluations.proof.linearisation_polynomial_eval,
        );

        // 5. Compute Openings using KZG10
        //
        // We merge the quotient polynomial using the `z_challenge` so the SRS
        // is linear in the circuit size `n`
        let quot = Self::compute_quotient_opening_poly(
            domain.size(),
            &t_1_poly,
            &t_2_poly,
            &t_3_poly,
            &t_4_poly,
            &z_challenge,
        );

        // Compute aggregate witness to polynomials evaluated at the evaluation
        // challenge `z`
        let aw_challenge: E::Fr =
            transcript.challenge_scalar(b"aggregate_witness");
        let aggregate_witness = Self::compute_aggregate_witness(
            &[
                quot,
                lin_poly,
                w_l_poly.clone(),
                w_r_poly.clone(),
                w_o_poly,
                w_4_poly.clone(),
                prover_key.permutation.left_sigma.0.clone(),
                prover_key.permutation.right_sigma.0.clone(),
                prover_key.permutation.out_sigma.0.clone(),
            ],
            &z_challenge,
            aw_challenge,
        );
        let w_z_comm = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key,
            &aggregate_witness,
            None,
            None,
        )?;

        // Compute aggregate witness to polynomials evaluated at the shifted
        // evaluation challenge
        let saw_challenge: E::Fr =
            transcript.challenge_scalar(b"aggregate_witness");
        let shifted_aggregate_witness = Self::compute_aggregate_witness(
            &[z_poly, w_l_poly, w_r_poly, w_4_poly],
            &(z_challenge * domain.element(1)),
            saw_challenge,
        );
        let w_zw_comm = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key,
            &shifted_aggregate_witness,
            None,
            None,
        )?;

        Ok(Proof {
            a_comm: w_l_poly_commit.0,
            b_comm: w_r_poly_commit.0,
            c_comm: w_o_poly_commit.0,
            d_comm: w_4_poly_commit.0,
            z_comm: z_poly_commit.0,
            t_1_comm: t_1_commit.0,
            t_2_comm: t_2_commit.0,
            t_3_comm: t_3_commit.0,
            t_4_comm: t_4_commit.0,
            w_z_comm: w_z_comm.0,
            w_zw_comm: w_zw_comm.0,
            evaluations: evaluations.proof,
        })
    }

    /// Proves a circuit is satisfied, then clears the witness variables
    /// If the circuit is not pre-processed, then the preprocessed circuit will
    /// also be computed.
    pub fn prove<E: PairingEngine<Fr = F>, P>(
        &mut self,
        commit_key: &Powers<E>,
    ) -> Result<Proof<E>, Error>
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
    {
        if self.prover_key.is_none() {
            // Preprocess circuit and store preprocessed circuit and transcript
            // in the Prover.
            self.prover_key = Some(self.cs.preprocess_prover(
                commit_key,
                &mut self.preprocessed_transcript,
            )?);
        }

        let prover_key = self.prover_key.as_ref().unwrap();

        let proof =
            self.prove_with_preprocessed::<E, P>(commit_key, prover_key)?;

        // Clear witness and reset composer variables
        self.clear_witness();

        Ok(proof)
    }
}

impl<F> Default for Prover<F>
where
    F: FftField,
{
    #[inline]
    fn default() -> Self {
        Prover::new(b"plonk")
    }
}
