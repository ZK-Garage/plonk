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
};
use ark_ec::{ModelParameters, TEModelParameters};
use ark_ff::{PrimeField, to_bytes};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
    UVPolynomial,
};
use core::marker::PhantomData;
use itertools::izip;
use ark_poly_commit::{QuerySet, PCRandomness};
use merlin::Transcript;

use ark_marlin::rng::FiatShamirRng;
use digest::Digest;
use blake2::Blake2s;

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

    /// Split `t(X)` poly into 4 n-sized polynomials.
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
    pub fn prove_with_preprocessed
    <D: Digest>(
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
        // let mut transcript = self.preprocessed_transcript.clone();

        pub const PROTOCOL_NAME: &[u8] = b"Plonk";
        let mut fs_rng = FiatShamirRng::<D>::from_seed(
            &to_bytes![
                &PROTOCOL_NAME
            ]
            .unwrap(),
        );

        // Append Public Inputs to the transcript
        // Add them in evaluations form since DensePolynomial doesn't implement to_bytes
        let pub_inputs = self.cs.get_pi().as_evals();
        fs_rng.absorb(&to_bytes![pub_inputs].unwrap());


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
        let (wire_commits, _) = PC::commit(commit_key, w_polys.iter(), None)
            .map_err(to_pc_error::<F, PC>)?;

        // Add witness polynomial commitments to transcript.
        fs_rng.absorb(&to_bytes![wire_commits].unwrap());

        // 2. Derive lookup polynomials

        // Generate table compression factor
        let zeta = F::rand(&mut fs_rng);
        fs_rng.absorb(&to_bytes![zeta].unwrap());

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
        fs_rng.absorb(&to_bytes![f_poly_commit].unwrap());


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
        fs_rng.absorb(&to_bytes![h_1_poly_commit, h_2_poly_commit].unwrap());


        // 3. Compute permutation polynomial
        //
        // Compute permutation challenge `beta`.
        let beta = F::rand(&mut fs_rng);
        fs_rng.absorb(&to_bytes![beta].unwrap());

        // Compute permutation challenge `gamma`.
        let gamma = F::rand(&mut fs_rng);
        fs_rng.absorb(&to_bytes![gamma].unwrap());

        // Compute permutation challenge `delta`.
        let delta = F::rand(&mut fs_rng);
        fs_rng.absorb(&to_bytes![delta].unwrap());

        let epsilon = F::rand(&mut fs_rng);
        fs_rng.absorb(&to_bytes![epsilon].unwrap());

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
        fs_rng.absorb(&to_bytes![z_poly_commit].unwrap());


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
        let pi_poly = self.cs.get_pi().into();

        // 4. Compute quotient polynomial
        // Compute quotient challenge; `alpha`, and gate-specific separation
        // challenges.

        let alpha = F::rand(&mut fs_rng);
        fs_rng.absorb(&to_bytes![alpha].unwrap());

        let range_sep_challenge = F::rand(&mut fs_rng);
        fs_rng.absorb(&to_bytes![range_sep_challenge].unwrap());

        let logic_sep_challenge = F::rand(&mut fs_rng);
        fs_rng.absorb(&to_bytes![logic_sep_challenge].unwrap());

        let fixed_base_sep_challenge = F::rand(&mut fs_rng);
        fs_rng.absorb(&to_bytes![fixed_base_sep_challenge].unwrap());

        let var_base_sep_challenge = F::rand(&mut fs_rng);
        fs_rng.absorb(&to_bytes![var_base_sep_challenge].unwrap());

        let lookup_sep_challenge = F::rand(&mut fs_rng);
        fs_rng.absorb(&to_bytes![lookup_sep_challenge].unwrap());

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
        fs_rng.absorb(&to_bytes![t_commits].unwrap());

        // 4. Compute linearisation polynomial
        // Compute evaluation challenge; `z`.

        let z_challenge = F::rand(&mut fs_rng);
        fs_rng.absorb(&to_bytes![z_challenge].unwrap());

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
            &t_1_poly,
            &t_2_poly,
            &t_3_poly,
            &t_4_poly,
            &z_poly,
            &z_2_poly,
            &f_poly,
            &h_1_poly,
            &h_2_poly,
            &table_poly,
        )?;

        // Add evaluations to transcript.
        // First wire evals

        fs_rng.absorb(&to_bytes![
            evaluations.wire_evals.a_eval, 
            evaluations.wire_evals.b_eval,
            evaluations.wire_evals.c_eval,
            evaluations.wire_evals.d_eval,
            evaluations.perm_evals.left_sigma_eval,
            evaluations.perm_evals.right_sigma_eval,
            evaluations.perm_evals.out_sigma_eval,
            evaluations.perm_evals.permutation_eval,
            evaluations.lookup_evals.f_eval,
            evaluations.lookup_evals.q_lookup_eval,
            evaluations.lookup_evals.z2_next_eval,
            evaluations.lookup_evals.h1_eval,
            evaluations.lookup_evals.h1_next_eval,
            evaluations.lookup_evals.h2_eval
        ].unwrap());

        // Third, all evals needed for custom gates
        evaluations
            .custom_evals
            .vals
            .iter()
            .for_each(|(_, eval)| {
                fs_rng.absorb(&to_bytes![eval].unwrap());
            });


        // 5. Compute Openings
        //
        let separation_challenge = F::rand(&mut fs_rng);

        let w_polys = [
            lin_poly.clone(),
            prover_key.permutation.left_sigma.0.clone(),
            prover_key.permutation.right_sigma.0.clone(),
            prover_key.permutation.out_sigma.0.clone(),
            f_poly,
            h_2_poly,
            table_poly.clone(),
            w_l_poly.clone(),
            w_r_poly.clone(),
            w_o_poly,
            w_4_poly.clone()
        ];

        // TODO: preprocess this commitments
        // adding PC to ProverKey introduces many changes
        // maybe there is a better place to store them or to introduce shared struct
        let (tmp_commits, _) = PC::commit(
            commit_key,
            &[
                label_polynomial!(lin_poly), // this can't be preprocessed but can be computed with MSM
                label_polynomial!(prover_key.permutation.left_sigma.0), 
                label_polynomial!(prover_key.permutation.right_sigma.0),
                label_polynomial!(prover_key.permutation.out_sigma.0),
                label_polynomial!(table_poly)

            ],
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;

        let w_commits = [
            tmp_commits[0].commitment().clone(),
            tmp_commits[1].commitment().clone(),
            tmp_commits[2].commitment().clone(),
            tmp_commits[3].commitment().clone(),
            f_poly_commit[0].commitment().clone(),
            h_2_poly_commit[0].commitment().clone(),
            tmp_commits[4].commitment().clone(),
            wire_commits[0].commitment().clone(),
            wire_commits[1].commitment().clone(),
            wire_commits[2].commitment().clone(),
            wire_commits[3].commitment().clone()
        ];

        let w_polys = w_polys.iter().enumerate().map(|(i, p)| {
            ark_poly_commit::LabeledPolynomial::new(format!("w_{}", i), p.clone(), None, None)
        }).collect::<Vec<ark_poly_commit::LabeledPolynomial<_, _>>>();

        let w_commits = w_commits.iter().enumerate().map(|(i, c)| {
            ark_poly_commit::LabeledCommitment::new(format!("w_{}", i), c.clone(), None)
        }).collect::<Vec<_>>();

        let saw_polys = [
            z_poly, 
            w_l_poly,
            w_r_poly,
            w_4_poly,
            h_1_poly,
            z_2_poly,
            table_poly
        ];

        let saw_commits = [
            z_poly_commit[0].commitment().clone(), 
            wire_commits[0].commitment().clone(),
            wire_commits[1].commitment().clone(),
            wire_commits[3].commitment().clone(),
            h_1_poly_commit[0].commitment().clone(),
            z_2_poly_commit[0].commitment().clone(),
            tmp_commits[4].commitment().clone()
        ];

        let saw_polys = saw_polys.iter().enumerate().map(|(i, p)| {
            ark_poly_commit::LabeledPolynomial::new(format!("saw_{}", i), p.clone(), None, None)
        }).collect::<Vec<_>>();

        let saw_commits = saw_commits.iter().enumerate().map(|(i, c)| {
            ark_poly_commit::LabeledCommitment::new(format!("saw_{}", i), c.clone(), None)
        }).collect::<Vec<_>>();

        let mut query_set = QuerySet::new();
        let z_label = String::from("z");
        let omega_z_label = String::from("omega_z");
        let omega_z_challenge = z_challenge * domain.element(1);

        for poly in &w_polys {
            query_set.insert((poly.label().clone(), (z_label.clone(), z_challenge)));
        }
        for poly in &saw_polys {
            query_set.insert((poly.label().clone(), (omega_z_label.clone(), omega_z_challenge)));
        }

        let rands = vec![PC::Randomness::empty(); w_commits.len() + saw_commits.len()];
        let batch_opening = PC::batch_open(
            commit_key,
            w_polys.iter().chain(saw_polys.iter()),
            w_commits.iter().chain(saw_commits.iter()),
            &query_set,
            separation_challenge,
            // w_rands.iter().chain(saw_rands.iter()),
            &rands,
            Some(&mut fs_rng),
        )
        .map_err(to_pc_error::<F, PC>)?;

        Ok(Proof {
            a_comm: wire_commits[0].commitment().clone(),
            b_comm: wire_commits[1].commitment().clone(),
            c_comm: wire_commits[2].commitment().clone(),
            d_comm: wire_commits[3].commitment().clone(),
            z_comm: z_poly_commit[0].commitment().clone(),
            f_comm: f_poly_commit[0].commitment().clone(),
            h_1_comm: h_1_poly_commit[0].commitment().clone(),
            h_2_comm: h_2_poly_commit[0].commitment().clone(),
            z_2_comm: z_2_poly_commit[0].commitment().clone(),
            t_1_comm: t_commits[0].commitment().clone(),
            t_2_comm: t_commits[1].commitment().clone(),
            t_3_comm: t_commits[2].commitment().clone(),
            t_4_comm: t_commits[3].commitment().clone(),
            evaluations,
            batch_opening
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

        let proof = self.prove_with_preprocessed::<Blake2s>(
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
