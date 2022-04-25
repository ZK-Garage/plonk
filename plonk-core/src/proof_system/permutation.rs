// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! PLONK Permutation Prover and Verifier Data

use crate::{
    error::Error,
    permutation::constants::{K1, K2, K3},
    proof_system::linearisation_poly::ProofEvaluations,
};
use ark_ff::FftField;
use ark_poly::{
    polynomial::univariate::DensePolynomial, EvaluationDomain, Evaluations,
    GeneralEvaluationDomain,
};
use ark_poly_commit::PCCommitment;
use ark_serialize::*;

/// Permutation Prover Key
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
    /// Left Permutation
    pub left_sigma: (DensePolynomial<F>, Evaluations<F>),

    /// Right Permutation
    pub right_sigma: (DensePolynomial<F>, Evaluations<F>),

    /// Output Permutation
    pub out_sigma: (DensePolynomial<F>, Evaluations<F>),

    /// Fourth Permutation
    pub fourth_sigma: (DensePolynomial<F>, Evaluations<F>),

    /// Linear Evaluations
    pub linear_evaluations: Evaluations<F>,
    /* Evaluations of f(x) = X
     * [XXX: Remove this and
     * benchmark if it makes a
     * considerable difference
     * -- These are just the
     * domain elements] */
}

impl<F> ProverKey<F>
where
    F: FftField,
{
    /// Computes permutation term of the quotient polynomial at the `i`th domain
    /// point.
    pub fn compute_quotient_i(
        &self,
        index: usize,
        w_l_i: F,
        w_r_i: F,
        w_o_i: F,
        w_4_i: F,
        z_i: F,
        z_i_next: F,
        alpha: F,
        l1_alpha_sq: F,
        beta: F,
        gamma: F,
    ) -> F {
        let a = self.compute_quotient_identity_range_check_i(
            index, w_l_i, w_r_i, w_o_i, w_4_i, z_i, alpha, beta, gamma,
        );
        let b = self.compute_quotient_copy_range_check_i(
            index, w_l_i, w_r_i, w_o_i, w_4_i, z_i_next, alpha, beta, gamma,
        );
        let c = self.compute_quotient_term_check_one_i(z_i, l1_alpha_sq);
        a + b + c
    }

    /// Computes the following:
    ///
    /// ```text
    /// (a(x) + beta * X + gamma) (b(X) + beta * k1 * X + gamma) (c(X) + beta *
    /// k2 * X + gamma)(d(X) + beta * k3 * X + gamma)z(X) * alpha
    /// ```
    fn compute_quotient_identity_range_check_i(
        &self,
        index: usize,
        w_l_i: F,
        w_r_i: F,
        w_o_i: F,
        w_4_i: F,
        z_i: F,
        alpha: F,
        beta: F,
        gamma: F,
    ) -> F {
        let x = self.linear_evaluations[index];
        (w_l_i + (beta * x) + gamma)
            * (w_r_i + (beta * K1::<F>() * x) + gamma)
            * (w_o_i + (beta * K2::<F>() * x) + gamma)
            * (w_4_i + (beta * K3::<F>() * x) + gamma)
            * z_i
            * alpha
    }

    /// Computes the following:
    ///
    /// ```text
    /// (a(x) + beta* Sigma1(X) + gamma) (b(X) + beta * Sigma2(X) + gamma) (c(X)
    /// + beta * Sigma3(X) + gamma)(d(X) + beta * Sigma4(X) + gamma) Z(X.omega) *
    /// alpha
    /// ```
    fn compute_quotient_copy_range_check_i(
        &self,
        index: usize,
        w_l_i: F,
        w_r_i: F,
        w_o_i: F,
        w_4_i: F,
        z_i_next: F,
        alpha: F,
        beta: F,
        gamma: F,
    ) -> F {
        let left_sigma_eval = self.left_sigma.1[index];
        let right_sigma_eval = self.right_sigma.1[index];
        let out_sigma_eval = self.out_sigma.1[index];
        let fourth_sigma_eval = self.fourth_sigma.1[index];
        let product = (w_l_i + (beta * left_sigma_eval) + gamma)
            * (w_r_i + (beta * right_sigma_eval) + gamma)
            * (w_o_i + (beta * out_sigma_eval) + gamma)
            * (w_4_i + (beta * fourth_sigma_eval) + gamma)
            * z_i_next
            * alpha;
        -product
    }

    /// Computes the following:
    ///
    /// ```text
    /// L_1(X)[Z(X) - 1]
    /// ```
    #[inline]
    fn compute_quotient_term_check_one_i(&self, z_i: F, l1_alpha_sq: F) -> F {
        (z_i - F::one()) * l1_alpha_sq
    }

    /// Computes the permutation term of the linearisation polynomial.
    pub fn compute_linearisation(
        &self,
        n: usize,
        z_challenge: F,
        (alpha, beta, gamma): (F, F, F),
        (a_eval, b_eval, c_eval, d_eval): (F, F, F, F),
        (sigma_1_eval, sigma_2_eval, sigma_3_eval): (F, F, F),
        z_eval: F,
        z_poly: &DensePolynomial<F>,
    ) -> Result<DensePolynomial<F>, Error> {
        let a = self.compute_lineariser_identity_range_check(
            (a_eval, b_eval, c_eval, d_eval),
            z_challenge,
            (alpha, beta, gamma),
            z_poly,
        );
        let b = self.compute_lineariser_copy_range_check(
            (a_eval, b_eval, c_eval),
            z_eval,
            sigma_1_eval,
            sigma_2_eval,
            sigma_3_eval,
            (alpha, beta, gamma),
            &self.fourth_sigma.0,
        );
        let domain = GeneralEvaluationDomain::new(n).ok_or(Error::InvalidEvalDomainSize {
            log_size_of_group: n.trailing_zeros(),
            adicity:
                <<F as FftField>::FftParams as ark_ff::FftParameters>::TWO_ADICITY,
        })?;
        let c = self.compute_lineariser_check_is_one(
            &domain,
            z_challenge,
            alpha.square(),
            z_poly,
        );
        Ok(&(&a + &b) + &c)
    }

    /// Computes the following:
    ///
    /// ```text
    /// (a_eval + beta * z_challenge + gamma)(b_eval + beta * K1 * z_challenge +
    /// gamma)(c_eval + beta * K2 * z_challenge + gamma) * alpha z(X)
    /// ```
    fn compute_lineariser_identity_range_check(
        &self,
        (a_eval, b_eval, c_eval, d_eval): (F, F, F, F),
        z_challenge: F,
        (alpha, beta, gamma): (F, F, F),
        z_poly: &DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        let beta_z = beta * z_challenge;

        // a_eval + beta * z_challenge + gamma
        let mut a_0 = a_eval + beta_z;
        a_0 += gamma;

        // b_eval + beta * K1 * z_challenge + gamma
        let beta_z_k1 = K1::<F>() * beta_z;
        let mut a_1 = b_eval + beta_z_k1;
        a_1 += gamma;

        // c_eval + beta * K2 * z_challenge + gamma
        let beta_z_k2 = K2::<F>() * beta_z;
        let mut a_2 = c_eval + beta_z_k2;
        a_2 += gamma;

        // d_eval + beta * K3 * z_challenge + gamma
        let beta_z_k3 = K3::<F>() * beta_z;
        let mut a_3 = d_eval + beta_z_k3;
        a_3 += gamma;

        let mut a = a_0 * a_1;
        a *= a_2;
        a *= a_3;
        a *= alpha; // (a_eval + beta * z_challenge + gamma)(b_eval + beta * K1 *
                    // z_challenge + gamma)(c_eval + beta * K2 * z_challenge + gamma)(d_eval
                    // + beta * K3 * z_challenge + gamma) * alpha
        z_poly * a // (a_eval + beta * z_challenge + gamma)(b_eval + beta * K1
                   // * z_challenge + gamma)(c_eval + beta * K2 * z_challenge +
                   // gamma) * alpha z(X)
    }

    /// Computes the following:
    ///
    /// ```text
    /// -(a_eval + beta * sigma_1 + gamma)(b_eval + beta * sigma_2 + gamma)
    /// (c_eval + beta * sigma_3 + gamma) * beta *z_eval * alpha^2 * Sigma_4(X)
    /// ```
    fn compute_lineariser_copy_range_check(
        &self,
        (a_eval, b_eval, c_eval): (F, F, F),
        z_eval: F,
        sigma_1_eval: F,
        sigma_2_eval: F,
        sigma_3_eval: F,
        (alpha, beta, gamma): (F, F, F),
        fourth_sigma_poly: &DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        // a_eval + beta * sigma_1 + gamma
        let beta_sigma_1 = beta * sigma_1_eval;
        let mut a_0 = a_eval + beta_sigma_1;
        a_0 += gamma;

        // b_eval + beta * sigma_2 + gamma
        let beta_sigma_2 = beta * sigma_2_eval;
        let mut a_1 = b_eval + beta_sigma_2;
        a_1 += gamma;

        // c_eval + beta * sigma_3 + gamma
        let beta_sigma_3 = beta * sigma_3_eval;
        let mut a_2 = c_eval + beta_sigma_3;
        a_2 += gamma;

        let beta_z_eval = beta * z_eval;

        let mut a = a_0 * a_1 * a_2;
        a *= beta_z_eval;
        a *= alpha; // (a_eval + beta * sigma_1 + gamma)(b_eval + beta * sigma_2 +
                    // gamma)(c_eval + beta * sigma_3 + gamma) * beta * z_eval * alpha

        fourth_sigma_poly * -a // -(a_eval + beta * sigma_1 + gamma)(b_eval +
                               // beta * sigma_2 + gamma) (c_eval + beta *
                               // sigma_3 + gamma) * beta * z_eval * alpha^2 *
                               // Sigma_4(X)
    }

    /// Computes the lineariser check.
    fn compute_lineariser_check_is_one(
        &self,
        domain: &GeneralEvaluationDomain<F>,
        z_challenge: F,
        alpha_sq: F,
        z_coeffs: &DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        let l_1_z = domain.evaluate_all_lagrange_coefficients(z_challenge)[0];
        z_coeffs * (l_1_z * alpha_sq)
    }
}

/// Permutation Verifier Key
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = "PCC: core::fmt::Debug"),
    Eq(bound = "PCC: Eq"),
    PartialEq(bound = "PCC: PartialEq")
)]
pub struct VerifierKey<PCC>
where
    PCC: PCCommitment + Default,
{
    /// Left Permutation Commitment
    pub left_sigma: PCC,

    /// Right Permutation Commitment
    pub right_sigma: PCC,

    /// Output Permutation Commitment
    pub out_sigma: PCC,

    /// Fourth Permutation Commitment
    pub fourth_sigma: PCC,
}

impl<PCC> VerifierKey<PCC>
where
    PCC: PCCommitment + Default,
{
    /// Computes the linearisation commitments.
    pub fn compute_linearisation_commitment<F: FftField>(
        &self,
        scalars: &mut Vec<F>,
        points: &mut Vec<PCC>,
        evaluations: &ProofEvaluations<F>,
        z_challenge: F,
        (alpha, beta, gamma): (F, F, F),
        l1_eval: F,
        z_comm: PCC,
    ) {
        let alpha_sq = alpha.square();

        // (a_eval + beta * z + gamma)(b_eval + beta * z * k1 +
        // gamma)(c_eval + beta * k2 * z + gamma)(d_eval + beta
        // * k3 * z + gamma) * alpha
        let x = {
            let beta_z = beta * z_challenge;
            let q_0 = evaluations.wire_evals.a_eval + beta_z + gamma;

            let beta_k1_z = beta * K1::<F>() * z_challenge;
            let q_1 = evaluations.wire_evals.b_eval + beta_k1_z + gamma;

            let beta_k2_z = beta * K2::<F>() * z_challenge;
            let q_2 = evaluations.wire_evals.c_eval + beta_k2_z + gamma;

            let beta_k3_z = beta * K3::<F>() * z_challenge;
            let q_3 =
                (evaluations.wire_evals.d_eval + beta_k3_z + gamma) * alpha;

            q_0 * q_1 * q_2 * q_3
        };

        // l1(z) * alpha^2
        let r = l1_eval * alpha_sq;

        scalars.push(x + r);
        points.push(z_comm);

        // -(a_eval + beta * sigma_1_eval + gamma)(b_eval + beta *
        // sigma_2_eval + gamma)(c_eval + beta * sigma_3_eval +
        // gamma) * alpha^2
        let y = {
            let beta_sigma_1 = beta * evaluations.perm_evals.left_sigma_eval;
            let q_0 = evaluations.wire_evals.a_eval + beta_sigma_1 + gamma;

            let beta_sigma_2 = beta * evaluations.perm_evals.right_sigma_eval;
            let q_1 = evaluations.wire_evals.b_eval + beta_sigma_2 + gamma;

            let beta_sigma_3 = beta * evaluations.perm_evals.out_sigma_eval;
            let q_2 = evaluations.wire_evals.c_eval + beta_sigma_3 + gamma;

            let q_3 = beta * evaluations.perm_evals.permutation_eval * alpha;

            -(q_0 * q_1 * q_2 * q_3)
        };

        scalars.push(y);
        points.push(self.fourth_sigma.clone());
    }
}
