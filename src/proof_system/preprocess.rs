// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Methods to preprocess the constraint system for use in a proof.

use crate::constraint_system::StandardComposer;
use crate::error::Error;
use crate::proof_system::{widget, ProverKey};
use ark_ec::{ModelParameters, PairingEngine};
use ark_ff::FftField;
use ark_ff::PrimeField;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain};
use ark_poly_commit::kzg10::{Powers, KZG10};
use ark_poly_commit::{LabeledPolynomial, PolynomialCommitment};
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
    q_arith: DensePolynomial<F>,
    q_range: DensePolynomial<F>,
    q_logic: DensePolynomial<F>,
    q_fixed_group_add: DensePolynomial<F>,
    q_variable_group_add: DensePolynomial<F>,
    left_sigma: DensePolynomial<F>,
    right_sigma: DensePolynomial<F>,
    out_sigma: DensePolynomial<F>,
    fourth_sigma: DensePolynomial<F>,
}

impl<F, P> StandardComposer<F, P>
where
    F: FftField,
    P: ModelParameters<BaseField = F>,
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
        self.q_arith.extend(zeroes_scalar.iter());
        self.q_range.extend(zeroes_scalar.iter());
        self.q_logic.extend(zeroes_scalar.iter());
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
            && self.q_arith.len() == k
            && self.q_range.len() == k
            && self.q_logic.len() == k
            && self.q_fixed_group_add.len() == k
            && self.q_variable_group_add.len() == k
            && self.w_l.len() == k
            && self.w_r.len() == k
            && self.w_o.len() == k
            && self.w_4.len() == k
        {
            Ok(())
        } else {
            Err(Error::MismatchedPolyLen)
        }
    }
}
impl<F, P> StandardComposer<F, P>
where
    F: FftField + PrimeField,
    P: ModelParameters<BaseField = F>,
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
        PC: PolynomialCommitment<F, DensePolynomial<F>>,
    {
        let (_, selectors, domain) =
            self.preprocess_shared(commit_key, transcript, _pc).unwrap();

        let domain_4n =
            GeneralEvaluationDomain::new(4 * domain.size()).unwrap();
        let q_m_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.q_m),
            domain_4n,
        );
        let q_l_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.q_l),
            domain_4n,
        );
        let q_r_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.q_r),
            domain_4n,
        );
        let q_o_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.q_o),
            domain_4n,
        );
        let q_c_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.q_c),
            domain_4n,
        );
        let q_4_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.q_4),
            domain_4n,
        );
        let q_arith_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.q_arith),
            domain_4n,
        );
        let q_range_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.q_range),
            domain_4n,
        );
        let q_logic_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.q_logic),
            domain_4n,
        );
        let q_fixed_group_add_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.q_fixed_group_add),
            domain_4n,
        );
        let q_variable_group_add_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.q_variable_group_add),
            domain_4n,
        );

        let left_sigma_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.left_sigma),
            domain_4n,
        );
        let right_sigma_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.right_sigma),
            domain_4n,
        );
        let out_sigma_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.out_sigma),
            domain_4n,
        );
        let fourth_sigma_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&selectors.fourth_sigma),
            domain_4n,
        );
        // XXX: Remove this and compute it on the fly
        let linear_eval_4n = Evaluations::from_vec_and_domain(
            domain_4n.coset_fft(&[F::zero(), F::one()]),
            domain_4n,
        );

        // Compute 4n evaluations for X^n -1
        let v_h_coset_4n =
            compute_vanishing_poly_over_coset(domain_4n, domain.size() as u64);

        Ok(ProverKey::from_polynomials_and_evals(
            domain.size(),
            (selectors.q_m, q_m_eval_4n),
            (selectors.q_l, q_l_eval_4n),
            (selectors.q_r, q_r_eval_4n),
            (selectors.q_o, q_o_eval_4n),
            (selectors.q_4, q_4_eval_4n),
            (selectors.q_c, q_c_eval_4n),
            (selectors.q_arith, q_arith_eval_4n),
            (selectors.q_range, q_range_eval_4n),
            (selectors.q_logic, q_logic_eval_4n),
            (selectors.q_fixed_group_add, q_fixed_group_add_eval_4n),
            (selectors.q_variable_group_add, q_variable_group_add_eval_4n),
            (selectors.left_sigma, left_sigma_eval_4n),
            (selectors.right_sigma, right_sigma_eval_4n),
            (selectors.out_sigma, out_sigma_eval_4n),
            (selectors.fourth_sigma, fourth_sigma_eval_4n),
            linear_eval_4n,
            v_h_coset_4n,
        ))
    }

    /// The verifier only requires the commitments in order to verify a
    /// [`Proof`](super::Proof) We can therefore speed up preprocessing for the
    /// verifier by skipping the FFTs needed to compute the 4n evaluations.
    pub fn preprocess_verifier<PC>(
        &mut self,
        commit_key: &PC::CommitterKey,
        transcript: &mut Transcript,
        _pc: PhantomData<PC>,
    ) -> Result<widget::VerifierKey<F, PC>, Error>
    where
        PC: PolynomialCommitment<F, DensePolynomial<F>>,
    {
        let (verifier_key, _, _) =
            self.preprocess_shared(commit_key, transcript, _pc).unwrap();
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
        ),
        Error,
    >
    where
        PC: PolynomialCommitment<F, DensePolynomial<F>>,
    {
        let domain = GeneralEvaluationDomain::new(self.circuit_size()).unwrap();

        // Check that the length of the wires is consistent.
        self.check_poly_same_len().unwrap();

        // 1. Pad circuit to a power of two
        self.pad(domain.size() as usize - self.n);

        let ifft = |q: &[F]| DensePolynomial {
            coeffs: domain.ifft(q),
        };

        let q_m_poly = ifft(&self.q_m);
        let q_r_poly = ifft(&self.q_r);
        let q_l_poly = ifft(&self.q_l);
        let q_o_poly = ifft(&self.q_o);
        let q_c_poly = ifft(&self.q_c);
        let q_4_poly = ifft(&self.q_4);
        let q_arith_poly = ifft(&self.q_arith);
        let q_range_poly = ifft(&self.q_range);
        let q_logic_poly = ifft(&self.q_logic);
        let q_fixed_group_add_poly = ifft(&self.q_fixed_group_add);
        let q_variable_group_add_poly = ifft(&self.q_variable_group_add);

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
                LabeledPolynomial::new(
                    "q_m_poly".to_owned(),
                    q_m_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "q_l_poly".to_owned(),
                    q_l_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "q_r_poly".to_owned(),
                    q_r_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "q_o_poly".to_owned(),
                    q_o_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "q_4_poly".to_owned(),
                    q_4_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "q_c_poly".to_owned(),
                    q_c_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "q_arith_poly".to_owned(),
                    q_arith_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "q_range_poly".to_owned(),
                    q_range_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "q_logic_poly".to_owned(),
                    q_logic_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "q_fixed_group_add_poly".to_owned(),
                    q_fixed_group_add_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "q_variable_group_add_poly".to_owned(),
                    q_variable_group_add_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "left_sigma_poly".to_owned(),
                    left_sigma_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "right_sigma_poly".to_owned(),
                    right_sigma_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "out_sigma_poly".to_owned(),
                    out_sigma_poly.clone(),
                    None,
                    None,
                ),
                LabeledPolynomial::new(
                    "fourth_sigma_poly".to_owned(),
                    fourth_sigma_poly.clone(),
                    None,
                    None,
                ),
            ]
            .iter(),
            None,
        )
        .unwrap();

        let verifier_key = widget::VerifierKey::from_polynomial_commitments(
            self.circuit_size(),
            commitments[0].commitment().clone(), // q_m_poly_commit.0,
            commitments[1].commitment().clone(), // q_l_poly_commit.0,
            commitments[2].commitment().clone(), // q_r_poly_commit.0,
            commitments[3].commitment().clone(), // q_o_poly_commit.0,
            commitments[4].commitment().clone(), // q_4_poly_commit.0,
            commitments[5].commitment().clone(), // q_c_poly_commit.0,
            commitments[6].commitment().clone(), // q_arith_poly_commit.0,
            commitments[7].commitment().clone(), // q_range_poly_commit.0,
            commitments[8].commitment().clone(), // q_logic_poly_commit.0,
            commitments[9].commitment().clone(), // q_fixed_group_add_poly_commit.0,
            commitments[10].commitment().clone(), // q_variable_group_add_poly_commit.0,
            commitments[11].commitment().clone(), // left_sigma_poly_commit.0,
            commitments[12].commitment().clone(), // right_sigma_poly_commit.0,
            commitments[13].commitment().clone(), // out_sigma_poly_commit.0,
            commitments[14].commitment().clone(), // fourth_sigma_poly_commit.0,
        );

        let selectors = SelectorPolynomials {
            q_m: q_m_poly.clone(),
            q_l: q_l_poly.clone(),
            q_r: q_r_poly.clone(),
            q_o: q_o_poly.clone(),
            q_c: q_c_poly.clone(),
            q_4: q_4_poly.clone(),
            q_arith: q_arith_poly.clone(),
            q_range: q_range_poly.clone(),
            q_logic: q_logic_poly.clone(),
            q_fixed_group_add: q_fixed_group_add_poly.clone(),
            q_variable_group_add: q_variable_group_add_poly.clone(),
            left_sigma: left_sigma_poly.clone(),
            right_sigma: right_sigma_poly.clone(),
            out_sigma: out_sigma_poly.clone(),
            fourth_sigma: fourth_sigma_poly.clone(),
        };

        // Add the circuit description to the transcript
        verifier_key.seed_transcript(transcript);

        Ok((verifier_key, selectors, domain))
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
    let coset_gen = F::multiplicative_generator().pow(&[poly_degree, 0, 0, 0]);
    let v_h: Vec<_> = (0..domain.size())
        .map(|i| {
            (coset_gen * group_gen.pow(&[poly_degree * i as u64, 0, 0, 0]))
                - F::one()
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
        F: FftField,
        P: ModelParameters<BaseField = F>,
    {
        let mut composer: StandardComposer<F, P> = StandardComposer::new();
        dummy_gadget(100, &mut composer);

        // Pad the circuit to next power of two
        let next_pow_2 = composer.n.next_power_of_two() as u64;
        composer.pad(next_pow_2 as usize - composer.n);

        let size = composer.n;
        assert!(size.is_power_of_two());
        assert!(composer.q_m.len() == size);
        assert!(composer.q_l.len() == size);
        assert!(composer.q_o.len() == size);
        assert!(composer.q_r.len() == size);
        assert!(composer.q_c.len() == size);
        assert!(composer.q_arith.len() == size);
        assert!(composer.q_range.len() == size);
        assert!(composer.q_logic.len() == size);
        assert!(composer.q_fixed_group_add.len() == size);
        assert!(composer.q_variable_group_add.len() == size);
        assert!(composer.w_l.len() == size);
        assert!(composer.w_r.len() == size);
        assert!(composer.w_o.len() == size);
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
