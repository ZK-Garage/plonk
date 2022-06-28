//! Parameters for a circuit

#![allow(missing_docs)]
use crate::{
    commitment::HomomorphicCommitment,
    constraint_system::{ecc::Point, StandardComposer},
    proof_system::GateConstraint,
};

use ark_ec::ModelParameters;
use ark_ff::PrimeField;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};

pub trait CircuitParameters: CanonicalSerialize + CanonicalDeserialize {
    type ScalarField: PrimeField;
    type EmbeddedCurve: EmbeddedCurve<Self>;
    type EmbeddedCurveParameters: ModelParameters;
    type PolynomialCommitment: HomomorphicCommitment<
        Self::ScalarField,
        Commitment = Self::Commitment,
        Proof = Self::Proof,
        VerifierKey = Self::VerifierKey,
        UniversalParams = Self::UniversalParams,
        CommitterKey = Self::CommitterKey,
    >;
    type Commitment: ark_poly_commit::PCCommitment + Default;
    type Proof: ark_poly_commit::PCProof;
    type VerifierKey: ark_poly_commit::PCVerifierKey;
    type UniversalParams;
    type CommitterKey;
    type FixedBaseScalarMul: GateConstraint<Self::ScalarField>;
    type CurveAddition: GateConstraint<Self::ScalarField>;
}

pub trait EmbeddedCurve<P: CircuitParameters>: Sized {
    fn identity(composer: &mut StandardComposer<P>) -> Point;
}

pub mod test {
    use super::CircuitParameters;
    use crate::proof_system::ecc::CurveAddition;
    use crate::proof_system::ecc::FixedBaseScalarMul;
    use crate::{
        commitment::{IPA, KZG10},
        proof_system::ecc::{SWEmbeddedCurve, TEEmbeddedCurve},
    };
    use ark_poly::univariate::DensePolynomial;
    use ark_poly_commit::PolynomialCommitment;
    use ark_serialize::{
        CanonicalDeserialize, CanonicalSerialize, Read, SerializationError,
        Write,
    };

    macro_rules! circuit_params_short_weierstrass {
        ($name:ident, $scalar:ty, $embeddedparameters:ty, $pc:ty) => {
            #[derive(CanonicalSerialize, CanonicalDeserialize)]
            pub struct $name;

            impl CircuitParameters for $name {
                type ScalarField = $scalar;
                type EmbeddedCurve = SWEmbeddedCurve<Self>;
                type EmbeddedCurveParameters = $embeddedparameters;
                type PolynomialCommitment = $pc;
                type Commitment =
                    <Self::PolynomialCommitment as PolynomialCommitment<
                        Self::ScalarField,
                        DensePolynomial<Self::ScalarField>,
                    >>::Commitment;
                type Proof =
                    <Self::PolynomialCommitment as PolynomialCommitment<
                        Self::ScalarField,
                        DensePolynomial<Self::ScalarField>,
                    >>::Proof;
                type VerifierKey =
                    <Self::PolynomialCommitment as PolynomialCommitment<
                        Self::ScalarField,
                        DensePolynomial<Self::ScalarField>,
                    >>::VerifierKey;
                type UniversalParams =
                    <Self::PolynomialCommitment as PolynomialCommitment<
                        Self::ScalarField,
                        DensePolynomial<Self::ScalarField>,
                    >>::UniversalParams;
                type CommitterKey =
                    <Self::PolynomialCommitment as PolynomialCommitment<
                        Self::ScalarField,
                        DensePolynomial<Self::ScalarField>,
                    >>::CommitterKey;
                type FixedBaseScalarMul =
                    FixedBaseScalarMul<Self, Self::EmbeddedCurve>;
                type CurveAddition = CurveAddition<Self, Self::EmbeddedCurve>;
            }
        };
    }

    macro_rules! circuit_params_twisted_edwards {
        ($name:ident, $scalar:ty, $embeddedparameters:ty, $pc:ty) => {
            #[derive(CanonicalSerialize, CanonicalDeserialize)]
            pub struct $name;

            impl CircuitParameters for $name {
                type ScalarField = $scalar;
                type EmbeddedCurve = TEEmbeddedCurve<Self>;
                type EmbeddedCurveParameters = $embeddedparameters;
                type PolynomialCommitment = $pc;
                type Commitment =
                    <Self::PolynomialCommitment as PolynomialCommitment<
                        Self::ScalarField,
                        DensePolynomial<Self::ScalarField>,
                    >>::Commitment;
                type Proof =
                    <Self::PolynomialCommitment as PolynomialCommitment<
                        Self::ScalarField,
                        DensePolynomial<Self::ScalarField>,
                    >>::Proof;
                type VerifierKey =
                    <Self::PolynomialCommitment as PolynomialCommitment<
                        Self::ScalarField,
                        DensePolynomial<Self::ScalarField>,
                    >>::VerifierKey;
                type UniversalParams =
                    <Self::PolynomialCommitment as PolynomialCommitment<
                        Self::ScalarField,
                        DensePolynomial<Self::ScalarField>,
                    >>::UniversalParams;
                type CommitterKey =
                    <Self::PolynomialCommitment as PolynomialCommitment<
                        Self::ScalarField,
                        DensePolynomial<Self::ScalarField>,
                    >>::CommitterKey;
                type FixedBaseScalarMul =
                    FixedBaseScalarMul<Self, Self::EmbeddedCurve>;
                type CurveAddition = CurveAddition<Self, Self::EmbeddedCurve>;
            }
        };
    }
    circuit_params_twisted_edwards!(
        Bls12_377_KZG,
        ark_bls12_377::Fr,
        ark_ed_on_bls12_377::EdwardsParameters,
        KZG10<ark_bls12_377::Bls12_377>
    );
    circuit_params_twisted_edwards!(Bls12_377_IPA, ark_bls12_377::Fr, ark_ed_on_bls12_377::EdwardsParameters, IPA<ark_bls12_377::G1Affine, blake2::Blake2s>);
    circuit_params_twisted_edwards!(
        Bls12_381_KZG,
        ark_bls12_381::Fr,
        ark_ed_on_bls12_381::EdwardsParameters,
        KZG10<ark_bls12_381::Bls12_381>
    );
    circuit_params_twisted_edwards!(Bls12_381_IPA, ark_bls12_381::Fr, ark_ed_on_bls12_381::EdwardsParameters, IPA<ark_bls12_381::G1Affine, blake2::Blake2s>);

    circuit_params_short_weierstrass!(
        BW6_761_KZG,
        ark_bw6_761::Fr,
        ark_bls12_377::g1::Parameters,
        KZG10<ark_bw6_761::BW6_761>
    );
    circuit_params_short_weierstrass!(Pallas_IPA, ark_pallas::Fr, ark_vesta::VestaParameters, IPA<ark_pallas::Affine, blake2::Blake2b>);
    circuit_params_short_weierstrass!(Vesta_IPA, ark_vesta::Fr, ark_pallas::PallasParameters, IPA<ark_vesta::Affine, blake2::Blake2b>);
}
