// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Prover-side of the PLONK Proving System

use crate::{
    commitment::HomomorphicCommitment,
    constraint_system::{StandardComposer, Variable},
    error::{to_pc_error, Error},
    label_polynomial,
    proof_system::{
        linearisation_poly, proof::Proof, quotient_poly, ProverKey,
    },
    transcript::TranscriptProtocol,
};
use ark_ec::{ModelParameters, TEModelParameters};
use ark_ff::PrimeField;
use ark_poly::{
    univariate::{DensePolynomial, SparsePolynomial},
    EvaluationDomain, GeneralEvaluationDomain, UVPolynomial,
};
use core::marker::PhantomData;
use core::ops::Add;
use num_traits::{One, Zero};
use crate::lookup::MultiSet;
use merlin::Transcript;
use rand::rngs::OsRng;

/// Abstraction structure designed to construct a circuit and generate
/// [`Proof`]s for it.
pub struct Prover<F, P, PC>
where
    F: PrimeField,
    P: ModelParameters<BaseField = F>,
    PC: HomomorphicCommitment<F>,
{
    /// Proving Key which is used to create proofs about a specific PLONK
    /// circuit.
    pub prover_key: Option<ProverKey<F>>,

    /// Circuit Description
    pub(crate) cs: StandardComposer<F, P>,

    /// Store the messages exchanged during the preprocessing stage.
    ///
    /// This is copied each time, we make a proof.
    pub preprocessed_transcript: Transcript,

    _phantom: PhantomData<PC>,
}
impl<F, P, PC> Prover<F, P, PC>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
    PC: HomomorphicCommitment<F>,
{
    /// Creates a new `Prover` instance.
    pub fn new(label: &'static [u8]) -> Self {
        Self {
            prover_key: None,
            cs: StandardComposer::new(),
            preprocessed_transcript: Transcript::new(label),
            _phantom: PhantomData::<PC>,
        }
    }

    /// Creates a new `Prover` object with some expected size.
    pub fn with_expected_size(label: &'static [u8], size: usize) -> Self {
        Self {
            prover_key: None,
            cs: StandardComposer::with_expected_size(size),
            preprocessed_transcript: Transcript::new(label),
            _phantom: PhantomData::<PC>,
        }
    }

    /// Returns a mutable copy of the underlying [`StandardComposer`].
    pub fn mut_cs(&mut self) -> &mut StandardComposer<F, P> {
        &mut self.cs
    }

    /// Returns the number of gates in the circuit thet the `Prover` actually
    /// stores inside.
    pub fn circuit_size(&self) -> usize {
        self.cs.circuit_size()
    }

    /// Preprocesses the underlying constraint system.
    pub fn preprocess(
        &mut self,
        commit_key: &PC::CommitterKey,
    ) -> Result<(), Error> {
        if self.prover_key.is_some() {
            return Err(Error::CircuitAlreadyPreprocessed);
        }
        let pk = self.cs.preprocess_prover(
            commit_key,
            &mut self.preprocessed_transcript,
            PhantomData::<PC>,
        )?;
        self.prover_key = Some(pk);
        Ok(())
    }

    /// Split `t(X)` poly into 4 polynomials.
    /// The first 3 polynomials have degree n, the 4th has degree n+6
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

    /// Adds to a given polynomial a blinder term of the form:
    /// (b0 + b1 X + ...+ bk X^k) Z_h
    /// where k is the hiding_degree and Z_h = X^n - 1, the vanishing
    /// polynomial.
    fn add_blinder(
        polynomial: &DensePolynomial<F>,
        n: usize,
        hiding_degree: usize,
    ) -> DensePolynomial<F> {
        if hiding_degree < n / 2 {
            let z_h: DensePolynomial<F> =
                SparsePolynomial::from_coefficients_slice(&[
                    (0, -F::one()),
                    (n, F::one()),
                ])
                .into();
            let rand_poly =
                DensePolynomial::from_coefficients_vec(vec![
                    F::rand(&mut OsRng);
                    hiding_degree + 1
                ]);
            let blinder_poly = &rand_poly * &z_h;
            polynomial + &blinder_poly
        } else {
            let mut sparse_blinder_vec =
                vec![(0, F::zero()); 2 * (hiding_degree + 1)];

            // Computes the multiplication of (b0 + b1X + ..+ bk X^k) (X^n -1)
            // = (- b0 -b1 X ... -bk X^k  ..., b0 X^n + b1 X^(n+1) + ... bk
            // X^(n+k) as long as k< n/2
            for i in 0..=hiding_degree {
                let random_blinder = F::rand(&mut OsRng);
                sparse_blinder_vec[i] = (i, -random_blinder);
                sparse_blinder_vec[hiding_degree + 1 + i] =
                    (n + i, random_blinder);
            }

            let blinder_poly =
                SparsePolynomial::from_coefficients_vec(sparse_blinder_vec);
            // panic!("The blinder poly is {:?}", blinder_poly);

            polynomial + &blinder_poly
        }
    }

    /// Creates a [`Proof]` that demonstrates that a circuit is satisfied.
    /// # Note
    /// If you intend to construct multiple [`Proof`]s with different witnesses,
    /// after calling this method, the user should then call
    /// [`Prover::clear_witness`].
    /// This is automatically done when [`Prover::prove`] is called.
    pub fn prove_with_preprocessed(
        &self,
        commit_key: &PC::CommitterKey,
        prover_key: &ProverKey<F>,
        _data: PhantomData<PC>,
    ) -> Result<Proof<F, PC>, Error> {
        let domain =
            GeneralEvaluationDomain::new(self.cs.circuit_size()).ok_or(Error::InvalidEvalDomainSize {
                log_size_of_group: self.cs.circuit_size().trailing_zeros(),
                adicity: <<F as ark_ff::FftField>::FftParams as ark_ff::FftParameters>::TWO_ADICITY,
            })?;
        let n = domain.size();

        // Since the caller is passing a pre-processed circuit
        // We assume that the Transcript has been seeded with the preprocessed
        // Commitments
        let mut transcript = self.preprocessed_transcript.clone();

        // 1. Compute witness Polynomials
        //
        // Convert Variables to scalars padding them to the
        // correct domain size.
        let pad = vec![F::zero(); n - self.cs.w_l.len()];
        let w_l_scalar = &[&self.to_scalars(&self.cs.w_l)[..], &pad].concat();
        let w_r_scalar = &[&self.to_scalars(&self.cs.w_r)[..], &pad].concat();
        let w_o_scalar = &[&self.to_scalars(&self.cs.w_o)[..], &pad].concat();
        let w_4_scalar = &[&self.to_scalars(&self.cs.w_4)[..], &pad].concat();

        // Ensure the lookup q_lookup selector is the correct size
        // for building the witness polynomial, f 
        let padded_q_lookup = [&self.cs.q_lookup[..], &pad].concat();

        // Witnesses are now in evaluation form, convert them to coefficients
        // so that we may commit to them.
        let mut w_l_poly =
            DensePolynomial::from_coefficients_vec(domain.ifft(w_l_scalar));
        let mut w_r_poly =
            DensePolynomial::from_coefficients_vec(domain.ifft(w_r_scalar));
        let mut w_o_poly =
            DensePolynomial::from_coefficients_vec(domain.ifft(w_o_scalar));
        let mut w_4_poly =
            DensePolynomial::from_coefficients_vec(domain.ifft(w_4_scalar));

        // Add blinders
        w_l_poly = Self::add_blinder(&w_l_poly, n, 1);
        w_r_poly = Self::add_blinder(&w_r_poly, n, 1);
        w_o_poly = Self::add_blinder(&w_o_poly, n, 1);
        w_4_poly = Self::add_blinder(&w_4_poly, n, 1);
        let w_polys = [
            label_polynomial!(w_l_poly),
            label_polynomial!(w_r_poly),
            label_polynomial!(w_o_poly),
            label_polynomial!(w_4_poly),
        ];

        // Commit to witness polynomials.
        let (w_commits, w_rands) = PC::commit(commit_key, w_polys.iter(), None)
            .map_err(to_pc_error::<F, PC>)?;

        // Add witness polynomial commitments to transcript.
        //transcript.append_commitments(&*w_commits, PhantomData::<PC>);
        transcript.append(b"w_l", w_commits[0].commitment());
        transcript.append(b"w_r", w_commits[1].commitment());
        transcript.append(b"w_o", w_commits[2].commitment());
        transcript.append(b"w_4", w_commits[3].commitment());
    
        // 2. Derive lookup polynomials
        
        // Generate table compression factor
        let zeta = transcript.challenge_scalar(b"zeta");
        transcript.append(b"zeta", &zeta);

        // Compress lookup table into vector of single elements
        let compressed_t_multiset = MultiSet::compress_four_arity(
            [
                &prover_key.lookup.table_1.0,
                &prover_key.lookup.table_2.0,
                &prover_key.lookup.table_3.0,
                &prover_key.lookup.table_4.0,
            ],
            zeta,
        );

        // Compute query table f
        // When q_lookup[i] is zero the wire value is replaced with a dummy
        // value Currently set as the first row of the public table
        // If q_lookup is one the wire values are preserved
        let f_1_scalar = w_l_scalar
            .iter()
            .zip(&padded_q_lookup)
            .map(|(w, s)| {
                *w * s + (F::one() - s) * compressed_t_multiset.0[0]
            })
            .collect::<Vec<F>>();
        let f_2_scalar = w_r_scalar
            .iter()
            .zip(&padded_q_lookup)
            .map(|(w, s)| *w * s)
            .collect::<Vec<F>>();
        let f_3_scalar = w_o_scalar
            .iter()
            .zip(&padded_q_lookup)
            .map(|(w, s)| *w * s)
            .collect::<Vec<F>>();
        let f_4_scalar = w_4_scalar
            .iter()
            .zip(&padded_q_lookup)
            .map(|(w, s)| *w * s)
            .collect::<Vec<F>>();

        // Compress all wires into a single vector
        let compressed_f_multiset = MultiSet::compress_four_arity(
            [
                &MultiSet::from(&f_1_scalar[..]),
                &MultiSet::from(&f_2_scalar[..]),
                &MultiSet::from(&f_3_scalar[..]),
                &MultiSet::from(&f_4_scalar[..]),
            ],
            zeta,
        );

        // Compute long query poly
        let f_poly = DensePolynomial::from_coefficients_vec(
            domain.ifft(&compressed_f_multiset.0),
        );

        // Commit to query polynomial
        let f_poly_commit = PC::commit(commit_key, &f_poly, None)?;

        // Add f_poly commitment to transcript
        transcript.append(b"f", &f_poly_commit.0);

        // Compute s, as the sorted and concatenated version of f and t
        let s = compressed_t_multiset
            .sorted_concat(&compressed_f_multiset)
            .unwrap();

        // Compute first and second halves of s, as h_1 and h_2
        let (h_1, h_2) = s.halve_alternating();

        // Compute h polys
        let h_1_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&h_1.0));
        let h_2_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&h_2.0));

        // Commit to h polys
        let h_1_poly_commit = PC::commit(commit_key, &h_1_poly, None).map_err(to_pc_error::<F, PC>);
        let h_2_poly_commit = PC::commit(commit_key, &h_2_poly, None).unwrap();

        // Add h polynomials to transcript
        transcript.append_commitment(b"h1", &h_1_poly_commit.0);
        transcript.append_commitment(b"h2", &h_2_poly_commit.0);

        // 3. Compute permutation polynomial 
        //
        // Compute permutation challenge `beta`.
        let beta = transcript.challenge_scalar(b"beta");
        transcript.append(b"beta", &beta);
        

        // Compute permutation challenge `gamma`.
        let gamma = transcript.challenge_scalar(b"gamma");
        transcript.append(b"gamma", &gamma);
        // Compute permutation challenge `delta`.
        let delta = transcript.challenge_scalar(b"delta");
        transcript.append(b"delta", &delta);

        // Compute permutation challenge `epsilon`.
        let epsilon = transcript.challenge_scalar(b"epsilon");
        transcript.append(b"epsilon", &epsilon);

        // Compute permutation challenge `theta`.
        let theta = transcript.challenge_scalar(b"theta");
        transcript.append(b"theta", &theta);
        
        // Challenges must be different
        assert!(beta != gamma, "challenges must be different");
        assert!(beta != delta, "challenges must be different");
        assert!(beta != epsilon, "challenges must be different");
        assert!(beta != theta, "challenges must be different");
        assert!(gamma != delta, "challenges must be different");
        assert!(gamma != epsilon, "challenges must be different");
        assert!(gamma != theta, "challenges must be different");
        assert!(delta != epsilon, "challenges must be different");
        assert!(delta != theta, "challenges must be different");
        assert!(epsilon != theta, "challenges must be different");

        let mut z_poly = DensePolynomial::from_coefficients_slice(
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

        // Add blinder for permutation poly
        z_poly = Self::add_blinder(&z_poly, n, 2);

        // Commit to permutation polynomial.
        let (z_poly_commit, _) =
            PC::commit(commit_key, &[label_polynomial!(z_poly)], None)
                .map_err(to_pc_error::<F, PC>)?;

        // Add permutation polynomial commitment to transcript.
        transcript.append(b"z", z_poly_commit[0].commitment());

        // Compute mega permutation polynomial.
        // Compute lookup permutation poly
        let p_poly = DensePolynomial::from_coefficients_slice(
            &self.cs.perm.compute_mega_permutation_poly(
                &domain,
                &compressed_f_multiset.0,
                &compressed_t_multiset.0,
                &h_1.0,
                &h_2.0,
                &delta,
                &epsilon,
            ),
        );

        // 3. Compute public inputs polynomial.
        let pi_poly = DensePolynomial::from_coefficients_vec(
            domain.ifft(&self.cs.construct_dense_pi_vec()),
        );

        // 4. Compute quotient polynomial
        //
        // Compute quotient challenge; `alpha`, and gate-specific separation
        // challenges.
        let alpha = transcript.challenge_scalar(b"alpha");
        transcript.append(b"alpha", &alpha);

        let range_sep_challenge =
            transcript.challenge_scalar(b"range separation challenge");
        transcript.append(b"range seperation challenge", &range_sep_challenge);

        let logic_sep_challenge =
            transcript.challenge_scalar(b"logic separation challenge");
        transcript.append(b"logic seperation challenge", &logic_sep_challenge);

        let fixed_base_sep_challenge =
            transcript.challenge_scalar(b"fixed base separation challenge");
        transcript.append(b"fixed base separation challenge", &fixed_base_sep_challenge);

        let var_base_sep_challenge =
            transcript.challenge_scalar(b"variable base separation challenge");
        transcript.append(b"variable base separation challenge", &var_base_sep_challenge);

        let lookup_sep_challenge =
            transcript.challenge_scalar(b"lookup separation challenge");
        transcript.append(b"lookup separation challenge", &lookup_sep_challenge);

        let t_poly = quotient_poly::compute::<F, P>(
            &domain,
            prover_key,
            &z_poly,
            &w_l_poly,
            &w_r_poly,
            &w_o_poly,
            &w_4_poly,
            &pi_poly,
            &f_poly,
            &alpha,
            &beta,
            &gamma,
            &zeta,
            &range_sep_challenge,
            &logic_sep_challenge,
            &fixed_base_sep_challenge,
            &var_base_sep_challenge,
            &lookup_sep_challenge,
        )?;

        // Split quotient polynomial into 4 degree `n` polynomials
        let (t_1_poly, t_2_poly, t_3_poly, t_4_poly) =
            self.split_tx_poly(n, &t_poly);

        // Commit to splitted quotient polynomial
        let (t_commits, _) = PC::commit(
            commit_key,
            &[
                label_polynomial!(t_1_poly),
                label_polynomial!(t_2_poly),
                label_polynomial!(t_3_poly),
                label_polynomial!(t_4_poly),
            ],
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;

        // Add quotient polynomial commitments to transcript
        transcript.append(b"t_1", t_commits[0].commitment());
        transcript.append(b"t_2", t_commits[1].commitment());
        transcript.append(b"t_3", t_commits[2].commitment());
        transcript.append(b"t_4", t_commits[3].commitment());

        // 4. Compute linearisation polynomial
        //
        // Compute evaluation challenge; `z`.
        let z_challenge = transcript.challenge_scalar(b"z");
        transcript.append_scalar(b"z", &z_challenge);

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
        )?;

        

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
            n,
            &t_1_poly,
            &t_2_poly,
            &t_3_poly,
            &t_4_poly,
            &z_challenge,
        );

        // Compute aggregate witness to polynomials evaluated at the evaluation
        // challenge `z`
        let aw_challenge: F = transcript.challenge_scalar(b"aggregate_witness");

        let aw_polys = [
            label_polynomial!(quot),
            label_polynomial!(lin_poly),
            label_polynomial!(prover_key.permutation.left_sigma.0.clone()),
            label_polynomial!(prover_key.permutation.right_sigma.0.clone()),
            label_polynomial!(prover_key.permutation.out_sigma.0.clone()),
        ];

        let (aw_commits, aw_rands) = PC::commit(commit_key, &aw_polys, None)
            .map_err(to_pc_error::<F, PC>)?;

        let aw_opening = PC::open(
            commit_key,
            aw_polys.iter().chain(w_polys.iter()),
            aw_commits.iter().chain(w_commits.iter()),
            &z_challenge,
            aw_challenge,
            aw_rands.iter().chain(w_rands.iter()),
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;

        let saw_challenge: F =
            transcript.challenge_scalar(b"aggregate_witness");

        let saw_polys = [
            label_polynomial!(z_poly),
            label_polynomial!(w_l_poly),
            label_polynomial!(w_r_poly),
            label_polynomial!(w_4_poly),
        ];

        let (saw_commits, saw_rands) = PC::commit(commit_key, &saw_polys, None)
            .map_err(to_pc_error::<F, PC>)?;

        let saw_opening = PC::open(
            commit_key,
            &saw_polys,
            &saw_commits,
            &(z_challenge * domain.element(1)),
            saw_challenge,
            &saw_rands,
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;

        Ok(Proof {
            a_comm: w_commits[0].commitment().clone(),
            b_comm: w_commits[1].commitment().clone(),
            c_comm: w_commits[2].commitment().clone(),
            d_comm: w_commits[3].commitment().clone(),
            z_comm: saw_commits[0].commitment().clone(),
            t_1_comm: t_commits[0].commitment().clone(),
            t_2_comm: t_commits[1].commitment().clone(),
            t_3_comm: t_commits[2].commitment().clone(),
            t_4_comm: t_commits[3].commitment().clone(),
            aw_opening,
            saw_opening,
            evaluations: evaluations.proof,
        })
    }

    /// Proves a circuit is satisfied, then clears the witness variables
    /// If the circuit is not pre-processed, then the preprocessed circuit will
    /// also be computed.
    pub fn prove(
        &mut self,
        commit_key: &PC::CommitterKey,
    ) -> Result<Proof<F, PC>, Error> {
        if self.prover_key.is_none() {
            // Preprocess circuit and store preprocessed circuit and transcript
            // in the Prover.
            self.prover_key = Some(self.cs.preprocess_prover(
                commit_key,
                &mut self.preprocessed_transcript,
                PhantomData::<PC>,
            )?);
        }

        let prover_key = self.prover_key.as_ref().unwrap();

        let proof = self.prove_with_preprocessed(
            commit_key,
            prover_key,
            PhantomData::<PC>,
        )?;

        // Clear witness and reset composer variables
        self.clear_witness();

        Ok(proof)
    }
}

impl<F, P, PC> Default for Prover<F, P, PC>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
    PC: HomomorphicCommitment<F>,
{
    #[inline]
    fn default() -> Self {
        Prover::new(b"plonk")
    }
}
