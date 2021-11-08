// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Methods to preprocess the constraint system for use in a proof

use crate::constraint_system::StandardComposer;
use crate::error::Error;
use crate::proof_system::{widget, ProverKey};
use crate::transcript::TranscriptWrapper;
use ark_ec::{PairingEngine, ProjectiveCurve, TEModelParameters};
use ark_ff::PrimeField;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{EvaluationDomain, Evaluations, GeneralEvaluationDomain};
use ark_poly_commit::kzg10::{Powers, KZG10};
use core::marker::PhantomData;
use num_traits::{One, Zero};

/// Struct that contains all of the selector and permutation [`Polynomial`]s in
/// PLONK.
pub(crate) struct SelectorPolynomials<F: PrimeField> {
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

impl<
        E: PairingEngine,
        T: ProjectiveCurve<BaseField = E::Fr>,
        P: TEModelParameters<BaseField = E::Fr>,
    > StandardComposer<E, T, P>
{
    /// Pads the circuit to the next power of two.
    ///
    /// # Note
    /// `diff` is the difference between circuit size and next power of two.
    fn pad(&mut self, diff: usize) {
        // Add a zero variable to circuit
        let zero_scalar = E::Fr::zero();
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

    /// These are the parts of preprocessing that the prover must compute
    /// Although the prover does not need the verification key, he must compute
    /// the commitments in order to seed the transcript, allowing both the
    /// prover and verifier to have the same view
    pub fn preprocess_prover(
        &mut self,
        commit_key: &Powers<E>,
        transcript: &mut TranscriptWrapper<E>,
    ) -> Result<ProverKey<E::Fr, P>, Error> {
        let (_, selectors, domain) =
            self.preprocess_shared(commit_key, transcript)?;

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
            domain_4n.coset_fft(&[E::Fr::zero(), E::Fr::one()]),
            domain_4n,
        );

        // Prover Key for arithmetic circuits
        let arithmetic_prover_key = widget::arithmetic::ProverKey {
            q_m: (selectors.q_m, q_m_eval_4n),
            q_l: (selectors.q_l.clone(), q_l_eval_4n.clone()),
            q_r: (selectors.q_r.clone(), q_r_eval_4n.clone()),
            q_o: (selectors.q_o, q_o_eval_4n),
            q_c: (selectors.q_c.clone(), q_c_eval_4n.clone()),
            q_4: (selectors.q_4, q_4_eval_4n),
            q_arith: (selectors.q_arith, q_arith_eval_4n),
        };

        // Prover Key for range circuits
        let range_prover_key = widget::range::ProverKey {
            q_range: (selectors.q_range, q_range_eval_4n),
        };

        // Prover Key for logic circuits
        let logic_prover_key = widget::logic::ProverKey {
            q_c: (selectors.q_c.clone(), q_c_eval_4n.clone()),
            q_logic: (selectors.q_logic, q_logic_eval_4n),
        };

        // Prover Key for ecc circuits
        let ecc_prover_key = widget::ecc::scalar_mul::fixed_base::ProverKey {
            q_l: (selectors.q_l, q_l_eval_4n),
            q_r: (selectors.q_r, q_r_eval_4n),
            q_c: (selectors.q_c, q_c_eval_4n),
            q_fixed_group_add: (
                selectors.q_fixed_group_add,
                q_fixed_group_add_eval_4n,
            ),
            _marker: PhantomData,
        };

        // Prover Key for permutation argument
        let permutation_prover_key = widget::permutation::ProverKey {
            left_sigma: (selectors.left_sigma, left_sigma_eval_4n),
            right_sigma: (selectors.right_sigma, right_sigma_eval_4n),
            out_sigma: (selectors.out_sigma, out_sigma_eval_4n),
            fourth_sigma: (selectors.fourth_sigma, fourth_sigma_eval_4n),
            linear_evaluations: linear_eval_4n,
        };

        // Prover Key for curve addition
        let curve_addition_prover_key =
            widget::ecc::curve_addition::ProverKey::new(
                selectors.q_variable_group_add,
                q_variable_group_add_eval_4n,
            );

        let prover_key = ProverKey {
            n: domain.size(),
            arithmetic: arithmetic_prover_key,
            logic: logic_prover_key,
            range: range_prover_key,
            permutation: permutation_prover_key,
            variable_base: curve_addition_prover_key,
            fixed_base: ecc_prover_key,
            // Compute 4n evaluations for X^n -1
            v_h_coset_4n: compute_vanishing_poly_over_coset(
                domain_4n,
                domain.size() as u64,
            ),
        };

        Ok(prover_key)
    }

    /// The verifier only requires the commitments in order to verify a
    /// [`Proof`](super::Proof) We can therefore speed up preprocessing for the
    /// verifier by skipping the FFTs needed to compute the 4n evaluations.
    pub fn preprocess_verifier(
        &mut self,
        commit_key: &Powers<E>,
        transcript: &mut TranscriptWrapper<E>,
    ) -> Result<widget::VerifierKey<E, P>, Error> {
        let (verifier_key, _, _) =
            self.preprocess_shared(commit_key, transcript)?;
        Ok(verifier_key)
    }

    /// Both the [`Prover`](super::Prover) and [`Verifier`](super::Verifier)
    /// must perform IFFTs on the selector polynomials and permutation
    /// polynomials in order to commit to them and have the same transcript
    /// view.
    fn preprocess_shared(
        &mut self,
        commit_key: &Powers<E>,
        transcript: &mut TranscriptWrapper<E>,
    ) -> Result<
        (
            widget::VerifierKey<E, P>,
            SelectorPolynomials<E::Fr>,
            GeneralEvaluationDomain<E::Fr>,
        ),
        Error,
    > {
        let domain = GeneralEvaluationDomain::new(self.circuit_size()).unwrap();

        // Check that the length of the wires is consistent.
        self.check_poly_same_len()?;

        // 1. Pad circuit to a power of two
        self.pad(domain.size() as usize - self.n);

        let q_m_poly: DensePolynomial<E::Fr> = DensePolynomial {
            coeffs: domain.ifft(&self.q_m),
        };
        let q_r_poly: DensePolynomial<E::Fr> = DensePolynomial {
            coeffs: domain.ifft(&self.q_r),
        };
        let q_l_poly: DensePolynomial<E::Fr> = DensePolynomial {
            coeffs: domain.ifft(&self.q_l),
        };
        let q_o_poly: DensePolynomial<E::Fr> = DensePolynomial {
            coeffs: domain.ifft(&self.q_o),
        };
        let q_c_poly: DensePolynomial<E::Fr> = DensePolynomial {
            coeffs: domain.ifft(&self.q_c),
        };
        let q_4_poly: DensePolynomial<E::Fr> = DensePolynomial {
            coeffs: domain.ifft(&self.q_4),
        };
        let q_arith_poly: DensePolynomial<E::Fr> = DensePolynomial {
            coeffs: domain.ifft(&self.q_arith),
        };
        let q_range_poly: DensePolynomial<E::Fr> = DensePolynomial {
            coeffs: domain.ifft(&self.q_range),
        };
        let q_logic_poly: DensePolynomial<E::Fr> = DensePolynomial {
            coeffs: domain.ifft(&self.q_logic),
        };
        let q_fixed_group_add_poly: DensePolynomial<E::Fr> = DensePolynomial {
            coeffs: domain.ifft(&self.q_fixed_group_add),
        };
        let q_variable_group_add_poly: DensePolynomial<E::Fr> =
            DensePolynomial {
                coeffs: domain.ifft(&self.q_variable_group_add),
            };

        // 2. Compute the sigma polynomials
        let (
            left_sigma_poly,
            right_sigma_poly,
            out_sigma_poly,
            fourth_sigma_poly,
        ) = self.perm.compute_sigma_polynomials(self.n, &domain);

        let q_m_poly_commit = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key, &q_m_poly, None, None,
        )?;

        let q_l_poly_commit = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key, &q_l_poly, None, None,
        )?;

        let q_r_poly_commit = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key, &q_r_poly, None, None,
        )?;

        let q_o_poly_commit = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key, &q_o_poly, None, None,
        )?;

        let q_c_poly_commit = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key, &q_c_poly, None, None,
        )?;

        let q_4_poly_commit = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key, &q_4_poly, None, None,
        )?;

        let q_arith_poly_commit = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key,
            &q_arith_poly,
            None,
            None,
        )?;

        let q_range_poly_commit = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key,
            &q_range_poly,
            None,
            None,
        )?;

        let q_logic_poly_commit = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key,
            &q_logic_poly,
            None,
            None,
        )?;

        let q_fixed_group_add_poly_commit =
            KZG10::<E, DensePolynomial<E::Fr>>::commit(
                commit_key,
                &q_fixed_group_add_poly,
                None,
                None,
            )?;

        let q_variable_group_add_poly_commit =
            KZG10::<E, DensePolynomial<E::Fr>>::commit(
                commit_key,
                &q_variable_group_add_poly,
                None,
                None,
            )?;

        let left_sigma_poly_commit =
            KZG10::<E, DensePolynomial<E::Fr>>::commit(
                commit_key,
                &left_sigma_poly,
                None,
                None,
            )?;

        let right_sigma_poly_commit =
            KZG10::<E, DensePolynomial<E::Fr>>::commit(
                commit_key,
                &right_sigma_poly,
                None,
                None,
            )?;

        let out_sigma_poly_commit = KZG10::<E, DensePolynomial<E::Fr>>::commit(
            commit_key,
            &out_sigma_poly,
            None,
            None,
        )?;

        let fourth_sigma_poly_commit =
            KZG10::<E, DensePolynomial<E::Fr>>::commit(
                commit_key,
                &fourth_sigma_poly,
                None,
                None,
            )?;

        let verifier_key: widget::VerifierKey<E, P> =
            widget::VerifierKey::from_polynomial_commitments(
                self.circuit_size(),
                q_m_poly_commit.0,
                q_l_poly_commit.0,
                q_r_poly_commit.0,
                q_o_poly_commit.0,
                q_4_poly_commit.0,
                q_c_poly_commit.0,
                q_arith_poly_commit.0,
                q_logic_poly_commit.0,
                q_range_poly_commit.0,
                q_fixed_group_add_poly_commit.0,
                q_variable_group_add_poly_commit.0,
                left_sigma_poly_commit.0,
                right_sigma_poly_commit.0,
                out_sigma_poly_commit.0,
                fourth_sigma_poly_commit.0,
            );

        let selectors = SelectorPolynomials {
            q_m: q_m_poly,
            q_l: q_l_poly,
            q_r: q_r_poly,
            q_o: q_o_poly,
            q_c: q_c_poly,
            q_4: q_4_poly,
            q_arith: q_arith_poly,
            q_range: q_range_poly,
            q_logic: q_logic_poly,
            q_fixed_group_add: q_fixed_group_add_poly,
            q_variable_group_add: q_variable_group_add_poly,
            left_sigma: left_sigma_poly,
            right_sigma: right_sigma_poly,
            out_sigma: out_sigma_poly,
            fourth_sigma: fourth_sigma_poly,
        };

        // Add the circuit description to the transcript
        verifier_key.seed_transcript(transcript);

        Ok((verifier_key, selectors, domain))
    }
}

/// Given that the domain size is `D`
/// This function computes the `D` evaluation points for
/// the vanishing polynomial of degree `n` over a coset
pub(crate) fn compute_vanishing_poly_over_coset<
    F: PrimeField,
    D: EvaluationDomain<F>,
>(
    domain: D,        // domain to evaluate over
    poly_degree: u64, // degree of the vanishing polynomial
) -> Evaluations<F, D> {
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
    use crate::constraint_system::helper::*;
    use ark_bls12_381::Bls12_381;
    use ark_ed_on_bls12_381::{
        EdwardsParameters as JubjubParameters,
        EdwardsProjective as JubjubProjective,
    };
    #[test]
    /// Tests that the circuit gets padded to the correct length
    /// XXX: We can do this test without dummy_gadget method
    fn test_pad() {
        let mut composer: StandardComposer<
            Bls12_381,
            JubjubProjective,
            JubjubParameters,
        > = StandardComposer::new();
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
}
