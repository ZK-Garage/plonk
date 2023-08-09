// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Prover-side of the PLONK Proving System

use crate::lookup::MultiSet;
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
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
    UVPolynomial,
};
use core::marker::PhantomData;
use itertools::izip;
use merlin::Transcript;

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

    /// Returns the smallest power of two needed for the curcuit.
    pub fn circuit_bound(&self) -> usize {
        self.cs.circuit_bound()
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

    /// Split `t(X)` poly into 8 n-sized polynomials.
    #[allow(clippy::type_complexity)] // NOTE: This is an ok type for internal use.
    fn split_tx_poly(
        &self,
        n: usize,
        t_x: &DensePolynomial<F>,
    ) -> [DensePolynomial<F>; 8] {
        let mut buf = t_x.coeffs.to_vec();
        buf.resize(n << 3, F::zero());

        [
            DensePolynomial::from_coefficients_vec(buf[0..n].to_vec()),
            DensePolynomial::from_coefficients_vec(buf[n..2 * n].to_vec()),
            DensePolynomial::from_coefficients_vec(buf[2 * n..3 * n].to_vec()),
            DensePolynomial::from_coefficients_vec(buf[3 * n..4 * n].to_vec()),
            DensePolynomial::from_coefficients_vec(buf[4 * n..5 * n].to_vec()),
            DensePolynomial::from_coefficients_vec(buf[5 * n..6 * n].to_vec()),
            DensePolynomial::from_coefficients_vec(buf[6 * n..7 * n].to_vec()),
            DensePolynomial::from_coefficients_vec(buf[7 * n..].to_vec()),
        ]
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

    /// Clears all data in the [`Prover`] instance.
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
            GeneralEvaluationDomain::new(self.cs.circuit_bound()).ok_or(Error::InvalidEvalDomainSize {
                log_size_of_group: self.cs.circuit_bound().trailing_zeros(),
                adicity: <<F as ark_ff::FftField>::FftParams as ark_ff::FftParameters>::TWO_ADICITY,
            })?;
        let n = domain.size();

        // Since the caller is passing a pre-processed circuit
        // We assume that the Transcript has been seeded with the preprocessed
        // Commitments
        let mut transcript = self.preprocessed_transcript.clone();

        // Append Public Inputs to the transcript
        transcript.append(b"pi", self.cs.get_pi());

        // 1. Compute witness Polynomials
        //
        // Convert Variables to scalars padding them to the
        // correct domain size.
        let pad = vec![F::zero(); n - self.cs.w_l.len()];
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
        transcript.append(b"w_l", w_commits[0].commitment());
        transcript.append(b"w_r", w_commits[1].commitment());
        transcript.append(b"w_o", w_commits[2].commitment());
        transcript.append(b"w_4", w_commits[3].commitment());

        // 2. Derive lookup polynomials

        // Generate table compression factor
        let zeta = transcript.challenge_scalar(b"zeta");
        transcript.append(b"zeta", &zeta);

        // Compress lookup table into vector of single elements
        let compressed_t_multiset = MultiSet::compress(
            &[
                prover_key.lookup.table_1.clone(),
                prover_key.lookup.table_2.clone(),
                prover_key.lookup.table_3.clone(),
                prover_key.lookup.table_4.clone(),
            ],
            zeta,
        );

        // Compute table poly
        let table_poly = DensePolynomial::from_coefficients_vec(
            domain.ifft(&compressed_t_multiset.0),
        );

        // Compute query table f
        // When q_lookup[i] is zero the wire value is replaced with a dummy
        //   value currently set as the first row of the public table
        // If q_lookup[i] is one the wire values are preserved
        // This ensures the ith element of the compressed query table
        //   is an element of the compressed lookup table even when
        //   q_lookup[i] is 0 so the lookup check will pass

        let q_lookup_pad = vec![F::zero(); n - self.cs.q_lookup.len()];
        let padded_q_lookup =
            &[self.cs.q_lookup.as_slice(), q_lookup_pad.as_slice()].concat();

        let mut f_scalars: Vec<MultiSet<F>> =
            vec![MultiSet::with_capacity(w_l_scalar.len()); 4];

        for (q_lookup, w_l, w_r, w_o, w_4) in izip!(
            padded_q_lookup,
            w_l_scalar,
            w_r_scalar,
            w_o_scalar,
            w_4_scalar,
        ) {
            if q_lookup.is_zero() {
                f_scalars[0].push(compressed_t_multiset.0[0]);
                f_scalars.iter_mut().skip(1).for_each(|f| f.push(F::zero()));
            } else {
                f_scalars[0].push(*w_l);
                f_scalars[1].push(*w_r);
                f_scalars[2].push(*w_o);
                f_scalars[3].push(*w_4);
            }
        }

        // Compress all wires into a single vector
        let compressed_f_multiset = MultiSet::compress(&f_scalars, zeta);

        // Compute query poly
        let f_poly = DensePolynomial::from_coefficients_vec(
            domain.ifft(&compressed_f_multiset.0),
        );

        // Add blinders to query polynomials
        // let f_poly = Self::add_blinder(&f_poly, n, 1);

        // Commit to query polynomial
        let (f_poly_commit, _) =
            PC::commit(commit_key, &[label_polynomial!(f_poly)], None)
                .map_err(to_pc_error::<F, PC>)?;

        // Add f_poly commitment to transcript
        transcript.append(b"f", f_poly_commit[0].commitment());

        // Compute s, as the sorted and concatenated version of f and t
        let (h_1, h_2) = compressed_t_multiset
            .combine_split(&compressed_f_multiset)
            .unwrap();

        // Compute h polys
        let h_1_poly =
            DensePolynomial::from_coefficients_vec(domain.ifft(&h_1.0));
        let h_2_poly =
            DensePolynomial::from_coefficients_vec(domain.ifft(&h_2.0));

        // Add blinders to h polynomials
        // let h_1_poly = Self::add_blinder(&h_1_poly, n, 1);
        // let h_2_poly = Self::add_blinder(&h_2_poly, n, 1);

        // Commit to h polys
        let (h_1_poly_commit, _) =
            PC::commit(commit_key, &[label_polynomial!(h_1_poly)], None)
                .map_err(to_pc_error::<F, PC>)?;
        let (h_2_poly_commit, _) =
            PC::commit(commit_key, &[label_polynomial!(h_2_poly)], None)
                .map_err(to_pc_error::<F, PC>)?;

        // Add h polynomials to transcript
        transcript.append(b"h1", h_1_poly_commit[0].commitment());
        transcript.append(b"h2", h_2_poly_commit[0].commitment());

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

        // Challenges must be different
        assert!(beta != gamma, "challenges must be different");
        assert!(beta != delta, "challenges must be different");
        assert!(beta != epsilon, "challenges must be different");
        assert!(gamma != delta, "challenges must be different");
        assert!(gamma != epsilon, "challenges must be different");
        assert!(delta != epsilon, "challenges must be different");

        let z_poly = self.cs.perm.compute_permutation_poly(
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
        );

        // Commit to permutation polynomial.
        let (z_poly_commit, _) =
            PC::commit(commit_key, &[label_polynomial!(z_poly)], None)
                .map_err(to_pc_error::<F, PC>)?;

        // Add permutation polynomial commitment to transcript.
        transcript.append(b"z", z_poly_commit[0].commitment());

        // Compute mega permutation polynomial.
        // Compute lookup permutation poly
        let z_2_poly = DensePolynomial::from_coefficients_slice(
            &self.cs.perm.compute_lookup_permutation_poly(
                &domain,
                &compressed_f_multiset.0,
                &compressed_t_multiset.0,
                &h_1.0,
                &h_2.0,
                delta,
                epsilon,
            ),
        );

        // TODO: Find strategy for blinding lookups
        // Add blinder for lookup permutation poly
        // z_2_poly = Self::add_blinder(&z_2_poly, n, 2);

        // Commit to lookup permutation polynomial.
        let (z_2_poly_commit, _) =
            PC::commit(commit_key, &[label_polynomial!(z_2_poly)], None)
                .map_err(to_pc_error::<F, PC>)?;

        // 3. Compute public inputs polynomial.
        let pi_poly = self.cs.get_pi().into_dense_poly(n);

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
        transcript.append(
            b"fixed base separation challenge",
            &fixed_base_sep_challenge,
        );

        let var_base_sep_challenge =
            transcript.challenge_scalar(b"variable base separation challenge");
        transcript.append(
            b"variable base separation challenge",
            &var_base_sep_challenge,
        );

        let lookup_sep_challenge =
            transcript.challenge_scalar(b"lookup separation challenge");
        transcript
            .append(b"lookup separation challenge", &lookup_sep_challenge);

        let t_poly = quotient_poly::compute::<F, P>(
            &domain,
            prover_key,
            &z_poly,
            &z_2_poly,
            &w_l_poly,
            &w_r_poly,
            &w_o_poly,
            &w_4_poly,
            &pi_poly,
            &f_poly,
            &table_poly,
            &h_1_poly,
            &h_2_poly,
            &alpha,
            &beta,
            &gamma,
            &delta,
            &epsilon,
            &zeta,
            &range_sep_challenge,
            &logic_sep_challenge,
            &fixed_base_sep_challenge,
            &var_base_sep_challenge,
            &lookup_sep_challenge,
        )?;

        let t_i_polys = self.split_tx_poly(n, &t_poly);
        // Commit to splitted quotient polynomial
        let (t_commits, _) = PC::commit(
            commit_key,
            &[
                label_polynomial!(t_i_polys[0]),
                label_polynomial!(t_i_polys[1]),
                label_polynomial!(t_i_polys[2]),
                label_polynomial!(t_i_polys[3]),
                label_polynomial!(t_i_polys[4]),
                label_polynomial!(t_i_polys[5]),
                label_polynomial!(t_i_polys[6]),
                label_polynomial!(t_i_polys[7]),
            ],
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;

        // Add quotient polynomial commitments to transcript
        transcript.append(b"t_1", t_commits[0].commitment());
        transcript.append(b"t_2", t_commits[1].commitment());
        transcript.append(b"t_3", t_commits[2].commitment());
        transcript.append(b"t_4", t_commits[3].commitment());
        transcript.append(b"t_5", t_commits[4].commitment());
        transcript.append(b"t_6", t_commits[5].commitment());
        transcript.append(b"t_7", t_commits[6].commitment());
        transcript.append(b"t_8", t_commits[7].commitment());

        // 4. Compute linearisation polynomial
        //
        // Compute evaluation challenge; `z`.
        let z_challenge = transcript.challenge_scalar(b"z");
        transcript.append(b"z", &z_challenge);

        let (lin_poly, evaluations) = linearisation_poly::compute::<F, P>(
            &domain,
            prover_key,
            &alpha,
            &beta,
            &gamma,
            &delta,
            &epsilon,
            &zeta,
            &range_sep_challenge,
            &logic_sep_challenge,
            &fixed_base_sep_challenge,
            &var_base_sep_challenge,
            &lookup_sep_challenge,
            &z_challenge,
            &w_l_poly,
            &w_r_poly,
            &w_o_poly,
            &w_4_poly,
            &t_i_polys[0],
            &t_i_polys[1],
            &t_i_polys[2],
            &t_i_polys[3],
            &t_i_polys[4],
            &t_i_polys[5],
            &t_i_polys[6],
            &t_i_polys[7],
            &z_poly,
            &z_2_poly,
            &f_poly,
            &h_1_poly,
            &h_2_poly,
            &table_poly,
        )?;

        // Add evaluations to transcript.
        // First wire evals
        transcript.append(b"a_eval", &evaluations.wire_evals.a_eval);
        transcript.append(b"b_eval", &evaluations.wire_evals.b_eval);
        transcript.append(b"c_eval", &evaluations.wire_evals.c_eval);
        transcript.append(b"d_eval", &evaluations.wire_evals.d_eval);

        // Second permutation evals
        transcript
            .append(b"left_sig_eval", &evaluations.perm_evals.left_sigma_eval);
        transcript.append(
            b"right_sig_eval",
            &evaluations.perm_evals.right_sigma_eval,
        );
        transcript
            .append(b"out_sig_eval", &evaluations.perm_evals.out_sigma_eval);
        transcript
            .append(b"perm_eval", &evaluations.perm_evals.permutation_eval);

        // Third lookup evals
        transcript.append(b"f_eval", &evaluations.lookup_evals.f_eval);
        transcript
            .append(b"q_lookup_eval", &evaluations.lookup_evals.q_lookup_eval);
        transcript.append(
            b"lookup_perm_eval",
            &evaluations.lookup_evals.z2_next_eval,
        );
        transcript.append(b"h_1_eval", &evaluations.lookup_evals.h1_eval);
        transcript
            .append(b"h_1_next_eval", &evaluations.lookup_evals.h1_next_eval);
        transcript.append(b"h_2_eval", &evaluations.lookup_evals.h2_eval);

        // Third, all evals needed for custom gates
        evaluations
            .custom_evals
            .vals
            .iter()
            .for_each(|(label, eval)| {
                let static_label = Box::leak(label.to_owned().into_boxed_str());
                transcript.append(static_label.as_bytes(), eval);
            });

        // 5. Compute Openings using KZG10
        //
        // We merge the quotient polynomial using the `z_challenge` so the SRS
        // is linear in the circuit size `n`

        // Compute aggregate witness to polynomials evaluated at the evaluation
        // challenge `z`
        let aw_challenge: F = transcript.challenge_scalar(b"aggregate_witness");

        // XXX: The quotient polynomials is used here and then in the
        // opening poly. It is being left in for now but it may not
        // be necessary. Warrants further investigation.
        // Ditto with the out_sigma poly.
        let aw_polys = [
            label_polynomial!(lin_poly),
            label_polynomial!(prover_key.permutation.left_sigma.0.clone()),
            label_polynomial!(prover_key.permutation.right_sigma.0.clone()),
            label_polynomial!(prover_key.permutation.out_sigma.0.clone()),
            label_polynomial!(f_poly),
            label_polynomial!(h_2_poly),
            label_polynomial!(table_poly),
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
            label_polynomial!(h_1_poly),
            label_polynomial!(z_2_poly),
            label_polynomial!(table_poly),
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
            f_comm: f_poly_commit[0].commitment().clone(),
            h_1_comm: h_1_poly_commit[0].commitment().clone(),
            h_2_comm: h_2_poly_commit[0].commitment().clone(),
            z_2_comm: z_2_poly_commit[0].commitment().clone(),
            t_1_comm: t_commits[0].commitment().clone(),
            t_2_comm: t_commits[1].commitment().clone(),
            t_3_comm: t_commits[2].commitment().clone(),
            t_4_comm: t_commits[3].commitment().clone(),
            t_5_comm: t_commits[4].commitment().clone(),
            t_6_comm: t_commits[5].commitment().clone(),
            t_7_comm: t_commits[6].commitment().clone(),
            t_8_comm: t_commits[7].commitment().clone(),
            aw_opening,
            saw_opening,
            evaluations,
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
