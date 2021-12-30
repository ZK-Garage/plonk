// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use ark_ec::{AffineCurve, ModelParameters, PairingEngine, TEModelParameters};
use ark_ff::{BigInteger, FftField, Field, FpParameters, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
    Polynomial, UVPolynomial,
};

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
    F: FftField,
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
pub(crate) fn from_embedded_curve_scalar<F, P>(
    embedded_scalar: <P as ModelParameters>::ScalarField,
) -> F
where
    F: FftField + PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    let scalar_repr = embedded_scalar.into_repr();
    let modulus = <<F as PrimeField>::Params as FpParameters>::MODULUS;
    if modulus.num_bits() >= scalar_repr.num_bits() {
        let s = <<F as PrimeField>::BigInt as BigInteger>::from_bits_le(
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
    F::from_le_bytes_mod_order(&scalar_repr.to_bytes_le())
}

/// Get a embedded curve scalar `P::ScalarField` from a scalar of the pariring
/// friendly curve. Panics if the pairing frindly curve scalar is greater than
/// the modulus of the embedded curve scalar field
#[allow(dead_code)]
pub(crate) fn to_embedded_curve_scalar<F, P>(pfc_scalar: F) -> P::ScalarField
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
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
        let m = <<F as PrimeField>::BigInt as BigInteger>::from_bits_le(
            &modulus.to_bits_le(),
        );
        assert!(scalar_repr < m,
            "The embedded scalar exceeds the capacity representation of the outter curve scalar");
    }
    P::ScalarField::from_le_bytes_mod_order(&scalar_repr.to_bytes_le())
}

use ark_ec::msm::VariableBaseMSM;
use ark_poly_commit::sonic_pc::SonicKZG10;
use ark_poly_commit::PolynomialCommitment;

pub trait HomomorphicCommitment<F>:
    PolynomialCommitment<F, DensePolynomial<F>>
where
    F: PrimeField,
{
    fn multi_scalar_mul(
        commitments: &[Self::Commitment],
        scalars: &[F],
    ) -> Self::Commitment;
}

type KZG10<E> = SonicKZG10<E, DensePolynomial<<E as PairingEngine>::Fr>>;
type KZG10Commitment<E> = <KZG10<E> as PolynomialCommitment<
    <E as PairingEngine>::Fr,
    DensePolynomial<<E as PairingEngine>::Fr>,
>>::Commitment;

impl<E> HomomorphicCommitment<E::Fr> for KZG10<E>
where
    E: PairingEngine,
{
    fn multi_scalar_mul(
        commitments: &[KZG10Commitment<E>],
        scalars: &[E::Fr],
    ) -> KZG10Commitment<E> {
        let scalars_repr = scalars
            .iter()
            .map(<E::Fr as PrimeField>::into_repr)
            .collect::<Vec<_>>();

        let points_repr = commitments.iter().map(|c| c.0).collect::<Vec<_>>();

        ark_poly_commit::kzg10::Commitment::<E>(
            VariableBaseMSM::multi_scalar_mul(&points_repr, &scalars_repr)
                .into(),
        )
    }
}

/// Computes a linear combination of the polynomial evaluations and polynomial
/// commitments provided a challenge.
// TODO: complete doc
pub fn linear_combination<F, H>(
    evals: &[F],
    commitments: &[H::Commitment],
    challenge: F,
) -> (H::Commitment, F)
where
    F: PrimeField,
    H: HomomorphicCommitment<F>,
{
    assert_eq!(evals.len(), commitments.len());
    let powers = powers_of(challenge).take(evals.len()).collect::<Vec<_>>();
    let combined_eval = evals
        .iter()
        .zip(powers.iter())
        .map(|(&eval, power)| eval * power)
        .sum();
    let combined_commitment = H::multi_scalar_mul(commitments, &powers);
    (combined_commitment, combined_eval)
}

/// Macro to quickly label polynomials
#[macro_export]
macro_rules! label_polynomial {
    ($poly:expr) => {
        ark_poly_commit::LabeledPolynomial::new(
            stringify!($poly).to_owned(),
            $poly.clone(),
            None,
            None,
        )
    };
}

/// Macro to quickly label polynomial commitments
#[macro_export]
macro_rules! label_commitment {
    ($comm:expr) => {
        ark_poly_commit::LabeledCommitment::new(
            stringify!($comm).to_owned(),
            $comm.clone(),
            None,
        )
    };
}
