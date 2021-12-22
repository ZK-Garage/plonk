// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-INFRA. All rights reserved.

use ark_ec::{AffineCurve, ModelParameters, PairingEngine, TEModelParameters};
use ark_ff::{BigInteger, FftField, Field, FpParameters, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
    Polynomial, UVPolynomial,
};
use ark_poly_commit::kzg10::Commitment;

/// Returns an iterator over increasing powers of the given `scalar` starting
/// at `0`.
#[inline]
pub fn powers_of<F>(scalar: F) -> impl Iterator<Item = F>
where
    F: Field,
{
    core::iter::successors(Some(F::one()), move |p| Some(*p * scalar))
}

/// Performs polynomial division by `(x - z)` with `x` indeterminant using
/// Ruffini's algorithm.
pub fn ruffini<F>(poly: DensePolynomial<F>, z: F) -> DensePolynomial<F>
where
    F: PrimeField,
{
    let mut quotient = Vec::with_capacity(poly.degree());
    let mut k = F::zero();

    // Reverse the results and use Ruffini's method to compute the quotient
    // The coefficients must be reversed as Ruffini's method
    // starts with the leading coefficient, while Polynomials
    // are stored in increasing order i.e. the leading coefficient is the
    // last element
    for coeff in poly.coeffs.into_iter().rev() {
        let t = coeff + k;
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

/// Evaluation Domain Extension Trait
pub trait EvaluationDomainExt<F>: EvaluationDomain<F>
where
    F: FftField,
{
    /// Returns the value of `log_2(self.size)`.
    fn log_size_of_group(&self) -> u32;

    /// Returns the inverse of the size in the field.
    fn size_inv(&self) -> F;

    /// Returns a fixed generator of the subgroup.
    fn group_gen(&self) -> F;

    /// Returns the inverse of the fixed generator of the subgroup.
    fn group_gen_inv(&self) -> F;

    /// Returns a fixed multiplicative generator of the finite field.
    fn generator_inv(&self) -> F;
}

impl<F> EvaluationDomainExt<F> for GeneralEvaluationDomain<F>
where
    F: FftField,
{
    #[inline]
    fn log_size_of_group(&self) -> u32 {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.log_size_of_group,
            GeneralEvaluationDomain::MixedRadix(domain) => {
                domain.log_size_of_group
            }
        }
    }

    #[inline]
    fn size_inv(&self) -> F {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.size_inv,
            GeneralEvaluationDomain::MixedRadix(domain) => domain.size_inv,
        }
    }

    #[inline]
    fn group_gen(&self) -> F {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.group_gen,
            GeneralEvaluationDomain::MixedRadix(domain) => domain.group_gen,
        }
    }

    #[inline]
    fn group_gen_inv(&self) -> F {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.group_gen_inv,
            GeneralEvaluationDomain::MixedRadix(domain) => domain.group_gen_inv,
        }
    }

    #[inline]
    fn generator_inv(&self) -> F {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.generator_inv,
            GeneralEvaluationDomain::MixedRadix(domain) => domain.generator_inv,
        }
    }
}

/// Get a pairing friendly curve scalar `E::Fr` from a scalar of the embedded
/// curve. Panics if the embedded scalar is greater than the modulus of the
/// pairing firendly curve scalar field
#[allow(dead_code)]
pub fn from_embedded_curve_scalar<E, P>(
    embedded_scalar: <P as ModelParameters>::ScalarField,
) -> E::Fr
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    let scalar_repr = embedded_scalar.into_repr();

    let modulus = <<E::Fr as PrimeField>::Params as FpParameters>::MODULUS;
    if modulus.num_bits() >= scalar_repr.num_bits() {
        let s = <<E::Fr as PrimeField>::BigInt as BigInteger>::from_bits_le(
            &scalar_repr.to_bits_le(),
        );
        assert!(s < modulus,
            "The embedded scalar exceeds the capacity representation of the outter curve scalar");
    } else {
        let m = <<P::ScalarField as PrimeField>::BigInt as BigInteger>::from_bits_le(
            &modulus.to_bits_le(),
        );
        assert!(scalar_repr < m,
            "The embedded scalar exceeds the capacity representation of the outter curve scalar");
    }
    E::Fr::from_le_bytes_mod_order(&scalar_repr.to_bytes_le())
}

/// Get a embedded curve scalar `P::ScalarField` from a scalar of the pariring
/// friendly curve. Panics if the pairing frindly curve scalar is greater than
/// the modulus of the embedded curve scalar field
#[allow(dead_code)]
pub(crate) fn to_embedded_curve_scalar<E, P>(
    pfc_scalar: E::Fr,
) -> P::ScalarField
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    let scalar_repr = pfc_scalar.into_repr();
    let modulus =
        <<P::ScalarField as PrimeField>::Params as FpParameters>::MODULUS;
    if modulus.num_bits() >= scalar_repr.num_bits() {
        let s = <<P::ScalarField as PrimeField>::BigInt as BigInteger>::from_bits_le(
            &scalar_repr.to_bits_le(),
        );
        assert!(s < modulus,
            "The embedded scalar exceeds the capacity representation of the outter curve scalar");
    } else {
        let m = <<E::Fr as PrimeField>::BigInt as BigInteger>::from_bits_le(
            &modulus.to_bits_le(),
        );
        assert!(scalar_repr < m,
            "The embedded scalar exceeds the capacity representation of the outter curve scalar");
    }
    P::ScalarField::from_le_bytes_mod_order(&scalar_repr.to_bytes_le())
}

/// Computes a linear combination of the polynomial evaluations and polynomial
/// commitments provided a challenge.
// TODO: complete doc
pub fn linear_combination<E>(
    evals: &[E::Fr],
    commitments: &[Commitment<E>],
    challenge: E::Fr,
) -> (Commitment<E>, E::Fr)
where
    E: PairingEngine,
{
    assert_eq!(evals.len(), commitments.len());
    let powers = powers_of(challenge).take(evals.len()).collect::<Vec<_>>();
    let combined_eval = evals
        .iter()
        .zip(powers.iter())
        .map(|(&eval, power)| eval * power)
        .sum();
    let combined_commitment = Commitment(
        commitments
            .iter()
            .zip(powers.iter())
            .map(|(commit, &power)| commit.0.mul(power))
            .sum::<E::G1Projective>()
            .into(),
    );
    (combined_commitment, combined_eval)
}
