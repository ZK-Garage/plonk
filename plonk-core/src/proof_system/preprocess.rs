// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Methods to preprocess the constraint system for use in a proof.

use crate::{
    commitment::HomomorphicCommitment,
    constraint_system::StandardComposer,
    error::{to_pc_error, Error},
    label_polynomial,
    lookup::PreprocessedLookupTable,
    proof_system::{widget, ProverKey},
};
use ark_ec::TEModelParameters;
use ark_ff::{FftField, PrimeField};
use ark_poly::{
    polynomial::univariate::DensePolynomial, EvaluationDomain, Evaluations,
    GeneralEvaluationDomain, UVPolynomial,
};
use core::marker::PhantomData;
use merlin::Transcript;

/// Struct that contains all of the selector and permutation [`Polynomial`]s in
/// PLONK.
///
/// [`Polynomial`]: DensePolynomial
pub struct SelectorPolynomials<F>
where
    F: FftField,
{
    q_m: DensePolynomial<F>,
    q_l: DensePolynomial<F>,
    q_r: DensePolynomial<F>,
    q_o: DensePolynomial<F>,
    q_c: DensePolynomial<F>,
    q_4: DensePolynomial<F>,
    q_hl: DensePolynomial<F>,
    q_hr: DensePolynomial<F>,
    q_h4: DensePolynomial<F>,
    q_arith: DensePolynomial<F>,
    q_range: DensePolynomial<F>,
    q_logic: DensePolynomial<F>,
    q_lookup: DensePolynomial<F>,
    q_fixed_group_add: DensePolynomial<F>,
    q_variable_group_add: DensePolynomial<F>,
    left_sigma: DensePolynomial<F>,
    right_sigma: DensePolynomial<F>,
    out_sigma: DensePolynomial<F>,
    fourth_sigma: DensePolynomial<F>,
}

impl<F, P> StandardComposer<F, P>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    /// Pads the circuit to the next power of two.
    ///
    /// # Note
    /// `diff` is the difference between circuit size and next power of two.
    fn pad(&mut self, diff: usize) {
        // Add a zero variable to circuit
        let zero_scalar = F::zero();
        let zero_var = self.zero_var();

        let zeroes_scalar = vec![zero_scalar; diff];
        let zeroes_var = vec![zero_var; diff];

        self.q_m.extend(zeroes_scalar.iter());
        self.q_l.extend(zeroes_scalar.iter());
        self.q_r.extend(zeroes_scalar.iter());
        self.q_o.extend(zeroes_scalar.iter());
        self.q_c.extend(zeroes_scalar.iter());
        self.q_4.extend(zeroes_scalar.iter());

        // add high degree selectors
        self.q_hl.extend(zeroes_scalar.iter());
        self.q_hr.extend(zeroes_scalar.iter());
        self.q_h4.extend(zeroes_scalar.iter());

        self.q_arith.extend(zeroes_scalar.iter());
        self.q_range.extend(zeroes_scalar.iter());
        self.q_logic.extend(zeroes_scalar.iter());
        self.q_lookup.extend(zeroes_scalar.iter());
        self.q_fixed_group_add.extend(zeroes_scalar.iter());
        self.q_variable_group_add.extend(zeroes_scalar.iter());

        self.w_l.extend(zeroes_var.iter());
        self.w_r.extend(zeroes_var.iter());
        self.w_o.extend(zeroes_var.iter());
        self.w_4.extend(zeroes_var.iter());

        self.n += diff;
    }

    /// Checks that all of the wires of the composer have the same
    /// length.
    fn check_poly_same_len(&self) -> Result<(), Error> {
        let k = self.q_m.len();

        if self.q_o.len() == k
            && self.q_l.len() == k
            && self.q_r.len() == k
            && self.q_c.len() == k
            && self.q_4.len() == k
            && self.q_hl.len() == k
            && self.q_hr.len() == k
            && self.q_h4.len() == k
            && self.q_arith.len() == k
            && self.q_range.len() == k
            && self.q_logic.len() == k
            && self.q_lookup.len() == k
            && self.q_fixed_group_add.len() == k
            && self.q_variable_group_add.len() == k
            && self.w_l.len() == k
            && self.w_r.len() == k
            && self.w_o.len() == k
            && self.w_4.len() == k
        {
            Ok(())
        } else {
            println!(
                "length: {:?}",
                [
                    self.q_l.len(),
                    self.q_r.len(),
                    self.q_c.len(),
                    self.q_4.len(),
                    self.q_hl.len(),
                    self.q_hr.len(),
                    self.q_h4.len(),
                    self.q_arith.len(),
                    self.q_range.len(),
                    self.q_logic.len(),
                    self.q_lookup.len(),
                    self.q_fixed_group_add.len(),
                    self.q_variable_group_add.len(),
                    self.w_l.len(),
                    self.w_r.len(),
                    self.w_o.len(),
                    self.w_4.len(),
                ]
                .as_ref()
            );

            Err(Error::MismatchedPolyLen)
        }
    }
}
impl<F, P> StandardComposer<F, P>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    /// These are the parts of preprocessing that the prover must compute
    /// Although the prover does not need the verification key, he must compute
    /// the commitments in order to seed the transcript, allowing both the
    /// prover and verifier to have the same view
    pub fn preprocess_prover<PC>(
        &mut self,
        commit_key: &PC::CommitterKey,
        transcript: &mut Transcript,
        _pc: PhantomData<PC>,
    ) -> Result<ProverKey<F>, Error>
    where
        PC: HomomorphicCommitment<F>,
    {
        let (_, selectors, domain, preprocessed_table) =
            self.preprocess_shared(commit_key, transcript, _pc)?;

        let domain_8n =
            GeneralEvaluationDomain::new(8 * domain.size()).ok_or(Error::InvalidEvalDomainSize {
                log_size_of_group: (8 * domain.size()).trailing_zeros(),
                adicity:
                    <<F as FftField>::FftParams as ark_ff::FftParameters>::TWO_ADICITY,
            })?;
        let q_m_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_m),
            domain_8n,
        );
        let q_l_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_l),
            domain_8n,
        );
        let q_r_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_r),
            domain_8n,
        );
        let q_o_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_o),
            domain_8n,
        );
        let q_c_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_c),
            domain_8n,
        );
        let q_4_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_4),
            domain_8n,
        );
        let q_hash_left_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_hl),
            domain_8n,
        );
        let q_hash_right_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_hr),
            domain_8n,
        );
        let q_hash_4_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_h4),
            domain_8n,
        );
        let q_arith_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_arith),
            domain_8n,
        );
        let q_range_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_range),
            domain_8n,
        );
        let q_logic_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_logic),
            domain_8n,
        );
        let q_lookup_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_lookup),
            domain_8n,
        );
        let q_fixed_group_add_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_fixed_group_add),
            domain_8n,
        );
        let q_variable_group_add_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.q_variable_group_add),
            domain_8n,
        );
        let left_sigma_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.left_sigma),
            domain_8n,
        );
        let right_sigma_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.right_sigma),
            domain_8n,
        );
        let out_sigma_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.out_sigma),
            domain_8n,
        );
        let fourth_sigma_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&selectors.fourth_sigma),
            domain_8n,
        );
        // XXX: Remove this and compute it on the fly
        let linear_eval_8n = Evaluations::from_vec_and_domain(
            domain_8n.coset_fft(&[F::zero(), F::one()]),
            domain_8n,
        );

        // Compute 8n evaluations for X^n -1
        let v_h_coset_8n =
            compute_vanishing_poly_over_coset(domain_8n, domain.size() as u64);

        Ok(ProverKey::from_polynomials_and_evals(
            domain.size(),
            (selectors.q_m, q_m_eval_8n),
            (selectors.q_l, q_l_eval_8n),
            (selectors.q_r, q_r_eval_8n),
            (selectors.q_o, q_o_eval_8n),
            (selectors.q_4, q_4_eval_8n),
            (selectors.q_c, q_c_eval_8n),
            (selectors.q_hl, q_hash_left_eval_8n),
            (selectors.q_hr, q_hash_right_eval_8n),
            (selectors.q_h4, q_hash_4_eval_8n),
            (selectors.q_arith, q_arith_eval_8n),
            (selectors.q_range, q_range_eval_8n),
            (selectors.q_logic, q_logic_eval_8n),
            (selectors.q_lookup, q_lookup_eval_8n),
            (selectors.q_fixed_group_add, q_fixed_group_add_eval_8n),
            (selectors.q_variable_group_add, q_variable_group_add_eval_8n),
            (selectors.left_sigma, left_sigma_eval_8n),
            (selectors.right_sigma, right_sigma_eval_8n),
            (selectors.out_sigma, out_sigma_eval_8n),
            (selectors.fourth_sigma, fourth_sigma_eval_8n),
            linear_eval_8n,
            v_h_coset_8n,
            preprocessed_table.t[0].0.clone(),
            preprocessed_table.t[1].0.clone(),
            preprocessed_table.t[2].0.clone(),
            preprocessed_table.t[3].0.clone(),
        ))
    }

    /// The verifier only requires the commitments in order to verify a
    /// [`Proof`](super::Proof) We can therefore speed up preprocessing for the
    /// verifier by skipping the FFTs needed to compute the 8n evaluations.
    pub fn preprocess_verifier<PC>(
        &mut self,
        commit_key: &PC::CommitterKey,
        transcript: &mut Transcript,
        _pc: PhantomData<PC>,
    ) -> Result<widget::VerifierKey<F, PC>, Error>
    where
        PC: HomomorphicCommitment<F>,
    {
        let (verifier_key, _, _, _) =
            self.preprocess_shared(commit_key, transcript, _pc)?;
        Ok(verifier_key)
    }

    /// Both the [`Prover`](super::Prover) and [`Verifier`](super::Verifier)
    /// must perform IFFTs on the selector polynomials and permutation
    /// polynomials in order to commit to them and have the same transcript
    /// view.
    #[allow(clippy::type_complexity)] // FIXME: Add struct for prover side (last two tuple items).
    fn preprocess_shared<PC>(
        &mut self,
        commit_key: &PC::CommitterKey,
        transcript: &mut Transcript,
        _pc: PhantomData<PC>,
    ) -> Result<
        (
            widget::VerifierKey<F, PC>,
            SelectorPolynomials<F>,
            GeneralEvaluationDomain<F>,
            PreprocessedLookupTable<F, PC>,
        ),
        Error,
    >
    where
        PC: HomomorphicCommitment<F>,
    {
        let domain = GeneralEvaluationDomain::new(self.circuit_bound()).ok_or(Error::InvalidEvalDomainSize {
            log_size_of_group: (self.circuit_bound()).trailing_zeros(),
            adicity:
                <<F as FftField>::FftParams as ark_ff::FftParameters>::TWO_ADICITY,
        })?;

        let preprocessed_table = PreprocessedLookupTable::<F, PC>::preprocess(
            &self.lookup_table,
            commit_key,
            domain.size() as u32,
        )
        .unwrap();

        // Check that the length of the wires is consistent.
        self.check_poly_same_len()?;

        // 1. Pad circuit to a power of two
        self.pad(domain.size() - self.n);

        let q_m_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_m));

        let q_r_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_r));

        let q_l_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_l));

        let q_o_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_o));

        let q_4_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_4));

        let q_c_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_c));

        let q_hl_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_hl));

        let q_hr_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_hr));

        let q_h4_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_h4));

        let q_arith_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_arith));

        let q_range_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_range));

        let q_logic_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_logic));

        let q_lookup_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(domain.ifft(&self.q_lookup));

        let q_fixed_group_add_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(
                domain.ifft(&self.q_fixed_group_add),
            );

        let q_variable_group_add_poly: DensePolynomial<F> =
            DensePolynomial::from_coefficients_vec(
                domain.ifft(&self.q_variable_group_add),
            );

        // 2. Compute the sigma polynomials
        let (
            left_sigma_poly,
            right_sigma_poly,
            out_sigma_poly,
            fourth_sigma_poly,
        ) = self.perm.compute_sigma_polynomials(self.n, &domain);

        let (commitments, _) = PC::commit(
            commit_key,
            [
                label_polynomial!(q_m_poly),
                label_polynomial!(q_l_poly),
                label_polynomial!(q_r_poly),
                label_polynomial!(q_o_poly),
                label_polynomial!(q_4_poly),
                label_polynomial!(q_c_poly),
                label_polynomial!(q_hl_poly),
                label_polynomial!(q_hr_poly),
                label_polynomial!(q_h4_poly),
                label_polynomial!(q_arith_poly),
                label_polynomial!(q_range_poly),
                label_polynomial!(q_logic_poly),
                label_polynomial!(q_lookup_poly),
                label_polynomial!(q_fixed_group_add_poly),
                label_polynomial!(q_variable_group_add_poly),
                label_polynomial!(left_sigma_poly),
                label_polynomial!(right_sigma_poly),
                label_polynomial!(out_sigma_poly),
                label_polynomial!(fourth_sigma_poly),
            ]
            .iter(),
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;

        let verifier_key = widget::VerifierKey::from_polynomial_commitments(
            self.n,
            commitments[0].commitment().clone(), // q_m
            commitments[1].commitment().clone(), // q_l
            commitments[2].commitment().clone(), // q_r
            commitments[3].commitment().clone(), // q_o
            commitments[4].commitment().clone(), // q_4
            commitments[5].commitment().clone(), // q_c
            commitments[6].commitment().clone(), // q_hl
            commitments[7].commitment().clone(), // q_hr
            commitments[8].commitment().clone(), // q_h4
            commitments[9].commitment().clone(), // q_arith
            commitments[10].commitment().clone(), // q_range
            commitments[11].commitment().clone(), // q_logic
            commitments[12].commitment().clone(), // q_lookup
            commitments[13].commitment().clone(), // q_fixed_group_add
            commitments[14].commitment().clone(), // q_variable_group_add
            commitments[15].commitment().clone(), // left_sigma
            commitments[16].commitment().clone(), // right_sigma
            commitments[17].commitment().clone(), // out_sigma
            commitments[18].commitment().clone(), // fourth_sigma
            preprocessed_table.t[0].1.clone(),
            preprocessed_table.t[1].1.clone(),
            preprocessed_table.t[2].1.clone(),
            preprocessed_table.t[3].1.clone(),
        );

        let selectors = SelectorPolynomials {
            q_m: q_m_poly,
            q_l: q_l_poly,
            q_r: q_r_poly,
            q_o: q_o_poly,
            q_c: q_c_poly,
            q_4: q_4_poly,
            q_hl: q_hl_poly,
            q_hr: q_hr_poly,
            q_h4: q_h4_poly,
            q_arith: q_arith_poly,
            q_range: q_range_poly,
            q_logic: q_logic_poly,
            q_lookup: q_lookup_poly,
            q_fixed_group_add: q_fixed_group_add_poly,
            q_variable_group_add: q_variable_group_add_poly,
            left_sigma: left_sigma_poly,
            right_sigma: right_sigma_poly,
            out_sigma: out_sigma_poly,
            fourth_sigma: fourth_sigma_poly,
        };

        // Add the circuit description to the transcript
        verifier_key.seed_transcript(transcript);

        Ok((verifier_key, selectors, domain, preprocessed_table))
    }
}

/// Given that the domain size is `D`
/// This function computes the `D` evaluation points for
/// the vanishing polynomial of degree `n` over a coset
pub fn compute_vanishing_poly_over_coset<F, D>(
    domain: D,        // domain to evaluate over
    poly_degree: u64, // degree of the vanishing polynomial
) -> Evaluations<F, D>
where
    F: FftField,
    D: EvaluationDomain<F>,
{
    assert!(
        (domain.size() as u64) > poly_degree,
        "domain_size = {}, poly_degree = {}",
        domain.size() as u64,
        poly_degree
    );
    let group_gen = domain.element(1);
    let coset_gen = F::multiplicative_generator().pow([poly_degree]);
    let v_h: Vec<_> = (0..domain.size())
        .map(|i| {
            (coset_gen * group_gen.pow([poly_degree * i as u64])) - F::one()
        })
        .collect();
    Evaluations::from_vec_and_domain(v_h, domain)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{batch_test_field_params, constraint_system::helper::*};
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;

    /// Tests that the circuit gets padded to the correct length.
    // FIXME: We can do this test without dummy_gadget method.
    fn test_pad<F, P>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
    {
        let mut composer: StandardComposer<F, P> = StandardComposer::new();
        dummy_gadget(100, &mut composer);

        // Pad the circuit to next power of two
        let next_pow_2 = composer.n.next_power_of_two() as u64;
        composer.pad(next_pow_2 as usize - composer.n);

        let size = composer.n;
        assert!(size.is_power_of_two());
        assert_eq!(composer.q_m.len(), size);
        assert_eq!(composer.q_l.len(), size);
        assert_eq!(composer.q_o.len(), size);
        assert_eq!(composer.q_r.len(), size);
        assert_eq!(composer.q_c.len(), size);
        assert_eq!(composer.q_hl.len(), size);
        assert_eq!(composer.q_hr.len(), size);
        assert_eq!(composer.q_h4.len(), size);
        assert_eq!(composer.q_arith.len(), size);
        assert_eq!(composer.q_range.len(), size);
        assert_eq!(composer.q_logic.len(), size);
        assert_eq!(composer.q_lookup.len(), size);
        assert_eq!(composer.q_fixed_group_add.len(), size);
        assert_eq!(composer.q_variable_group_add.len(), size);
        assert_eq!(composer.w_l.len(), size);
        assert_eq!(composer.w_r.len(), size);
        assert_eq!(composer.w_o.len(), size);
    }

    // Bls12-381 tests
    batch_test_field_params!(
        [test_pad],
        [] => (
            Bls12_381,
            ark_ed_on_bls12_381::EdwardsParameters
        )
    );

    // Bls12-377 tests
    batch_test_field_params!(
        [test_pad],
        [] => (
            Bls12_377,
            ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
