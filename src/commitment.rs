//! Useful commitment stuff
use ark_ec::{msm::VariableBaseMSM, PairingEngine};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::sonic_pc::SonicKZG10;
use ark_poly_commit::PolynomialCommitment;

/// A homomorphic polynomial commitment
pub trait HomomorphicCommitment<F>:
    PolynomialCommitment<F, DensePolynomial<F>>
where
    F: PrimeField,
{
    /// Combine a linear combination of homomorphic commitments
    fn multi_scalar_mul(
        commitments: &[Self::Commitment],
        scalars: &[F],
    ) -> Self::Commitment;
}

/// The Default KZG-style commitment scheme
pub type KZG10<E> = SonicKZG10<E, DensePolynomial<<E as PairingEngine>::Fr>>;
/// A single KZG10 commitment
pub type KZG10Commitment<E> = <KZG10<E> as PolynomialCommitment<
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
    let powers = crate::util::powers_of(challenge)
        .take(evals.len())
        .collect::<Vec<_>>();
    let combined_eval = evals
        .iter()
        .zip(powers.iter())
        .map(|(&eval, power)| eval * power)
        .sum();
    let combined_commitment = H::multi_scalar_mul(commitments, &powers);
    (combined_commitment, combined_eval)
}
