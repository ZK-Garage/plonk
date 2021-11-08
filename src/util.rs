// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use ark_ec::{AffineCurve, ModelParameters, PairingEngine, TEModelParameters};
use ark_ff::{BigInteger, FftField, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, GeneralEvaluationDomain, Polynomial,
    UVPolynomial,
};
use ark_poly_commit::kzg10::Commitment;

/// Returns a vector of scalars of increasing powers of x from x^0 to x^d.
pub(crate) fn powers_of<F: PrimeField>(
    scalar: &F,
    max_degree: usize,
) -> Vec<F> {
    let mut powers = Vec::with_capacity(max_degree + 1);
    powers.push(F::one());
    for i in 1..=max_degree {
        powers.push(powers[i - 1] * scalar);
    }
    powers
}

pub fn ruffini<F: PrimeField>(
    poly: DensePolynomial<F>,
    z: F,
) -> DensePolynomial<F> {
    let mut quotient: Vec<F> = Vec::with_capacity(poly.degree());
    let mut k = F::zero();

    // Reverse the results and use Ruffini's method to compute the quotient
    // The coefficients must be reversed as Ruffini's method
    // starts with the leading coefficient, while Polynomials
    // are stored in increasing order i.e. the leading coefficient is the
    // last element
    for coeff in poly.coeffs.iter().rev() {
        let t = *coeff + k;
        quotient.push(t);
        k = z * t;
    }

    // Pop off the last element, it is the remainder term
    // For PLONK, we only care about perfect factors
    quotient.pop();

    // Reverse the results for storage in the Polynomial struct
    quotient.reverse();
    DensePolynomial::from_coefficients_vec(quotient)
}

pub fn get_domain_attrs<F: FftField>(
    domain: &GeneralEvaluationDomain<F>,
    attr: &'static str,
) -> F {
    match domain {
        GeneralEvaluationDomain::MixedRadix(domain) => match attr {
            "group_gen" => domain.group_gen,
            "group_gen_inv" => domain.group_gen_inv,
            "generator_inv" => domain.generator_inv,
            "size_as_fe" => domain.size_as_field_element,
            "size_inv" => domain.size_inv,
            _ => unreachable!(),
        },
        GeneralEvaluationDomain::Radix2(domain) => match attr {
            "group_gen" => domain.group_gen,
            "group_gen_inv" => domain.group_gen_inv,
            "generator_inv" => domain.generator_inv,
            "size_as_fe" => domain.size_as_field_element,
            "size_inv" => domain.size_inv,
            _ => unreachable!(),
        },
        _ => unreachable!(),
    }
}

// Get a representation of an embedded curve scalar as a scalar of the pairing friendly curve
pub fn from_embedded_curve_scalar<
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
>(
    embedded_scalar: <P as ModelParameters>::ScalarField,
) -> E::Fr {
    let scalar_repr = embedded_scalar.into_repr();
    E::Fr::from_le_bytes_mod_order(&scalar_repr.to_bytes_le())
}

// Computes a linear combination of the polynomial evaluations and polynomial
// commitments provided a challenge.
// TODO complete doc
pub fn linear_combination<E: PairingEngine>(
    evals: &[E::Fr],
    commitments: &[Commitment<E>],
    challenge: E::Fr,
) -> (Commitment<E>, E::Fr) {
    assert_eq!(evals.len(), commitments.len());
    // Generate a challenge to generate a linear combination of the proofs
    let powers: Vec<E::Fr> =
        crate::util::powers_of(&challenge, evals.len() - 1);
    let combined_eval: E::Fr = evals
        .iter()
        .zip(powers.iter())
        .map(|(&eval, power)| eval * power)
        .sum();

    let combined_commitment: Commitment<E> = Commitment(
        commitments
            .iter()
            .zip(powers.iter())
            .map(|(commit, &power)| commit.0.mul(power))
            .sum::<E::G1Projective>()
            .into(),
    );
    (combined_commitment, combined_eval)
}
