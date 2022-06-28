// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Tools & traits for PLONK circuits

use crate::{
    error::{to_pc_error, Error},
    parameters::CircuitParameters,
    prelude::StandardComposer,
    proof_system::{
        pi::PublicInputs, Proof, Prover, ProverKey, Verifier, VerifierKey,
    },
};
use ark_poly_commit::PolynomialCommitment;
use ark_serialize::*;

/// Collection of structs/objects that the Verifier will use in order to
/// de/serialize data needed for Circuit proof verification.
/// This structure can be seen as a link between the [`Circuit`] public input
/// positions and the [`VerifierKey`] that the Verifier needs to use.
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = "VerifierKey<P>: core::fmt::Debug"),
    Eq(bound = "VerifierKey<P>: Eq"),
    PartialEq(bound = "VerifierKey<P>: PartialEq")
)]
pub struct VerifierData<P>
where
    P: CircuitParameters,
{
    /// Verifier Key
    pub key: VerifierKey<P>,
    /// Public Input
    pub pi: PublicInputs<P::ScalarField>,
}

impl<P> VerifierData<P>
where
    P: CircuitParameters,
{
    /// Creates a new `VerifierData` from a [`VerifierKey`] and the public
    /// input positions of the circuit that it represents.
    pub fn new(key: VerifierKey<P>, pi: PublicInputs<P::ScalarField>) -> Self {
        Self { key, pi }
    }

    /// Returns a reference to the contained [`VerifierKey`].
    pub fn key(&self) -> &VerifierKey<P> {
        &self.key
    }

    /// Returns a reference to the contained Public Input .
    pub fn pi(&self) -> &PublicInputs<P::ScalarField> {
        &self.pi
    }
}

/// Trait that should be implemented for any circuit function to provide to it
/// the capabilities of automatically being able to generate, and verify proofs
/// as well as compile the circuit.
///
/// # Example
///
/// ```rust,no_run
/// use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
/// use ark_ec::PairingEngine;
/// use ark_ec::models::twisted_edwards_extended::GroupAffine;
/// use ark_ec::{TEModelParameters, AffineCurve, ProjectiveCurve};
/// use ark_ed_on_bls12_381::{
///     EdwardsAffine as JubJubAffine, EdwardsParameters as JubJubParameters,
///     EdwardsProjective as JubJubProjective, Fr as JubJubScalar,
/// };
/// use ark_ff::{FftField, PrimeField, BigInteger, ToConstraintField};
/// use plonk_core::circuit::{Circuit, verify_proof};
/// use plonk_core::constraint_system::StandardComposer;
/// use plonk_core::error::{to_pc_error,Error};
/// use ark_poly::polynomial::univariate::DensePolynomial;
/// use ark_poly_commit::{PolynomialCommitment, sonic_pc::SonicKZG10};
/// use plonk_core::prelude::*;
/// use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
/// use num_traits::{Zero, One};
/// use rand_core::OsRng;
///
/// fn main() -> Result<(), Error> {
/// // Implements a circuit that checks:
/// // 1) a + b = c where C is a PI
/// // 2) a <= 2^6
/// // 3) b <= 2^5
/// // 4) a * b = d where D is a PI
/// // 5) JubJub::GENERATOR * e(JubJubScalar) = f where F is a PI
/// #[derive(derivative::Derivative)]
///    #[derivative(Debug(bound = ""), Default(bound = ""))]
/// pub struct TestCircuit<
/// F: FftField + PrimeField,
/// P: TEModelParameters<BaseField = F>,
/// > {
/// a: F,
/// b: F,
/// c: F,
/// d: F,
/// e: P::ScalarField,
/// f: GroupAffine<P>,
/// }

/// impl<P, EmbeddedBaseField, EmbeddedCurveParameters> Circuit<P>
/// for TestCircuit<EmbeddedBaseField, EmbeddedCurveParameters>
/// where
/// EmbeddedBaseField: PrimeField,
/// EmbeddedCurveParameters:
///     TEModelParameters<BaseField = EmbeddedBaseField>,
/// P: CircuitParameters<
///     ScalarField = EmbeddedBaseField,
///     EmbeddedCurve = TEEmbeddedCurve<P>,
///     EmbeddedCurveParameters = EmbeddedCurveParameters
/// >,
///    {
///        const CIRCUIT_ID: [u8; 32] = [0xff; 32];
///
///        fn gadget(
///            &mut self,
///            composer: &mut StandardComposer<P>,
///        ) -> Result<(), Error> {
///            let a = composer.add_input(self.a);
///            let b = composer.add_input(self.b);
///            let zero = composer.zero_var();
///
///            // Make first constraint a + b = c (as public input)
///            composer.arithmetic_gate(|gate| {
///                gate.witness(a, b, Some(zero))
///                    .add(P::ScalarField::one(), P::ScalarField::one())
///                    .pi(-self.c)
///            });
///
///            // Check that a and b are in range
///            composer.range_gate(a, 1 << 6);
///            composer.range_gate(b, 1 << 5);
///            // Make second constraint a * b = d
///            composer.arithmetic_gate(|gate| {
///                gate.witness(a, b,
/// Some(zero)).mul(P::ScalarField::one()).pi(-self.d)            });
///            let e = composer
///                .add_input(from_embedded_curve_scalar::<EmbeddedBaseField,
/// EmbeddedCurveParameters>(self.e));             let (x, y) =
/// EmbeddedCurveParameters::AFFINE_GENERATOR_COEFFS;             let generator
/// = GroupAffine::new(x, y);             let scalar_mul_result =
///                composer.fixed_base_scalar_mul(e, generator);
///            // Apply the constrain
///            composer.assert_equal_public_point(scalar_mul_result, self.f);
///            Ok(())
///        }
///
///        fn padded_circuit_size(&self) -> usize {
///            1 << 11
///        }
///    }
///
/// // Generate CRS
/// type PC = SonicKZG10::<Bls12_381,DensePolynomial<BlsScalar>>;
/// type P = plonk_core::parameters::test::Bls12_381_KZG;
/// let pp = PC::setup(
///     1 << 10, None, &mut OsRng
///  )?;
///
/// let mut circuit = TestCircuit::<BlsScalar, JubJubParameters>::default();
/// // Compile the circuit
/// let (pk_p, vk) = circuit.compile(&pp)?;
///
/// let (x, y) = JubJubParameters::AFFINE_GENERATOR_COEFFS;
/// let generator: GroupAffine<JubJubParameters> = GroupAffine::new(x, y);
/// let point_f_pi: GroupAffine<JubJubParameters> = AffineCurve::mul(
///     &generator,
///     JubJubScalar::from(2u64).into_repr(),
/// )
/// .into_affine();
/// // Prover POV
/// let (proof, pi) = {
///     let mut circuit: TestCircuit<BlsScalar, JubJubParameters> = TestCircuit
/// {         a: BlsScalar::from(20u64),
///         b: BlsScalar::from(5u64),
///         c: BlsScalar::from(25u64),
///         d: BlsScalar::from(100u64),
///         e: JubJubScalar::from(2u64),
///         f: point_f_pi,
///     };
///     circuit.gen_proof(&pp, pk_p, b"Test")
/// }?;
///
/// let verifier_data = VerifierData::new(vk, pi);
/// // Test serialisation for verifier_data
/// let mut verifier_data_bytes = Vec::new();
/// verifier_data.serialize(&mut verifier_data_bytes).unwrap();
///
/// let deserialized_verifier_data: VerifierData<P> =
///     VerifierData::deserialize(verifier_data_bytes.as_slice()).unwrap();
///
/// assert!(deserialized_verifier_data == verifier_data);
///
/// // Verifier POV
/// verify_proof::<P>(
///     &pp,
///     verifier_data.key,
///     &proof,
///     &verifier_data.pi,
///     b"Test",
/// )
/// }
/// ```
pub trait Circuit<P>
where
    P: CircuitParameters,
{
    /// Circuit identifier associated constant.
    const CIRCUIT_ID: [u8; 32];

    /// Gadget implementation used to fill the composer.
    fn gadget(
        &mut self,
        composer: &mut StandardComposer<P>,
    ) -> Result<(), Error>;

    /// Compiles the circuit by using a function that returns a `Result`
    /// with the [`ProverKey`], [`VerifierKey`] and the circuit size.
    #[allow(clippy::type_complexity)] // NOTE: Clippy is too harsh here.
    fn compile(
        &mut self,
        u_params: &P::UniversalParams,
    ) -> Result<(ProverKey<P::ScalarField>, VerifierKey<P>), Error> {
        // Setup PublicParams
        let circuit_size = self.padded_circuit_size();
        let (ck, _) =
            P::PolynomialCommitment::trim(u_params, circuit_size, 0, None)
                .map_err(to_pc_error::<P>)?;

        //Generate & save `ProverKey` with some random values.
        let mut prover = Prover::<P>::new(b"CircuitCompilation");
        self.gadget(prover.mut_cs())?;
        prover.preprocess(&ck)?;

        // Generate & save `VerifierKey` with some random values.
        let mut verifier = Verifier::new(b"CircuitCompilation");
        self.gadget(verifier.mut_cs())?;
        verifier.preprocess(&ck)?;
        Ok((
            prover
                .prover_key
                .expect("Unexpected error. Missing ProverKey in compilation"),
            verifier
                .verifier_key
                .expect("Unexpected error. Missing VerifierKey in compilation"),
        ))
    }

    /// Generates a proof using the provided [`ProverKey`] and
    /// [`ark_poly_commit::PCUniversalParams`]. Returns a
    /// [`crate::proof_system::Proof`] and the [`PublicInputs`].
    fn gen_proof(
        &mut self,
        u_params: &P::UniversalParams,
        prover_key: ProverKey<P::ScalarField>,
        transcript_init: &'static [u8],
    ) -> Result<(Proof<P>, PublicInputs<P::ScalarField>), Error> {
        let circuit_size = self.padded_circuit_size();
        let (ck, _) =
            P::PolynomialCommitment::trim(u_params, circuit_size, 0, None)
                .map_err(to_pc_error::<P>)?;
        // New Prover instance
        let mut prover = Prover::new(transcript_init);
        // Fill witnesses for Prover
        self.gadget(prover.mut_cs())?;
        // Add ProverKey to Prover
        prover.prover_key = Some(prover_key);
        let pi = prover.cs.get_pi().clone();

        Ok((prover.prove(&ck)?, pi))
    }

    /// Returns the Circuit size padded to the next power of two.
    fn padded_circuit_size(&self) -> usize;
}

/// Verifies a proof using the provided `CircuitInputs` & `VerifierKey`
/// instances.
pub fn verify_proof<P>(
    u_params: &P::UniversalParams,
    plonk_verifier_key: VerifierKey<P>,
    proof: &Proof<P>,
    public_inputs: &PublicInputs<P::ScalarField>,
    transcript_init: &'static [u8],
) -> Result<(), Error>
where
    P: CircuitParameters,
{
    let mut verifier: Verifier<P> = Verifier::new(transcript_init);
    let padded_circuit_size = plonk_verifier_key.padded_circuit_size();
    verifier.verifier_key = Some(plonk_verifier_key);
    let (_, vk) =
        P::PolynomialCommitment::trim(u_params, padded_circuit_size, 0, None)
            .map_err(to_pc_error::<P>)?;

    verifier.verify(proof, &vk, public_inputs)
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::parameters::test::*;
    use crate::{
        constraint_system::StandardComposer,
        proof_system::ecc::TEEmbeddedCurve, util,
    };
    use ark_ec::{
        models::TEModelParameters, twisted_edwards_extended::GroupAffine,
        AffineCurve, ProjectiveCurve,
    };
    use ark_ff::{FftField, PrimeField};
    use ark_poly_commit::PolynomialCommitment;
    use rand_core::OsRng;

    // Implements a circuit that checks:
    // 1) a + b = c where C is a PI
    // 2) a <= 2^6
    // 3) b <= 2^5
    // 4) a * b = d where D is a PI
    // 5) JubJub::GENERATOR * e(JubJubScalar) = f where F is a PI
    #[derive(derivative::Derivative)]
    #[derivative(Debug(bound = ""), Default(bound = ""))]
    pub struct TestCircuit<
        F: FftField + PrimeField,
        P: TEModelParameters<BaseField = F>,
    > {
        a: F,
        b: F,
        c: F,
        d: F,
        e: P::ScalarField,
        f: GroupAffine<P>,
    }

    impl<P, EmbeddedBaseField, EmbeddedCurveParameters> Circuit<P>
        for TestCircuit<EmbeddedBaseField, EmbeddedCurveParameters>
    where
        EmbeddedBaseField: PrimeField,
        EmbeddedCurveParameters:
            TEModelParameters<BaseField = EmbeddedBaseField>,
        P: CircuitParameters<
            ScalarField = EmbeddedBaseField,
            EmbeddedCurve = TEEmbeddedCurve<P>,
            EmbeddedCurveParameters = EmbeddedCurveParameters,
        >,
    {
        const CIRCUIT_ID: [u8; 32] = [0xff; 32];

        fn gadget(
            &mut self,
            composer: &mut StandardComposer<P>,
        ) -> Result<(), Error> {
            let a = composer.add_input(self.a);
            let b = composer.add_input(self.b);
            let zero = composer.zero_var;

            // Make first constraint a + b = c (as public input)
            composer.arithmetic_gate(|gate| {
                gate.witness(a, b, Some(zero))
                    .add(P::ScalarField::one(), P::ScalarField::one())
                    .pi(-self.c)
            });

            // Check that a and b are in range
            composer.range_gate(a, 1 << 6);
            composer.range_gate(b, 1 << 5);
            // Make second constraint a * b = d
            composer.arithmetic_gate(|gate| {
                gate.witness(a, b, Some(zero))
                    .mul(P::ScalarField::one())
                    .pi(-self.d)
            });
            let e = composer.add_input(util::from_embedded_curve_scalar::<
                P::ScalarField,
                EmbeddedCurveParameters,
            >(self.e));
            let (x, y) = EmbeddedCurveParameters::AFFINE_GENERATOR_COEFFS;
            let generator = GroupAffine::new(x, y);
            let scalar_mul_result =
                composer.fixed_base_scalar_mul(e, generator);

            // Apply the constrain
            composer.assert_equal_public_point(scalar_mul_result, self.f);
            Ok(())
        }

        fn padded_circuit_size(&self) -> usize {
            1 << 9
        }
    }

    fn test_full<P, EmbeddedBaseField, EmbeddedCurveParameters>(
    ) -> Result<(), Error>
    where
        EmbeddedBaseField: PrimeField,
        EmbeddedCurveParameters:
            TEModelParameters<BaseField = EmbeddedBaseField>,
        P: CircuitParameters<
            ScalarField = EmbeddedBaseField,
            EmbeddedCurve = TEEmbeddedCurve<P>,
            EmbeddedCurveParameters = EmbeddedCurveParameters,
        >,
        VerifierData<P>: PartialEq,
    {
        // Generate CRS
        let pp = P::PolynomialCommitment::setup(1 << 10, None, &mut OsRng)
            .map_err(to_pc_error::<P>)?;

        let mut circuit =
            TestCircuit::<P::ScalarField, EmbeddedCurveParameters>::default();

        // Compile the circuit
        let (pk, vk) = circuit.compile(&pp)?;

        let (x, y) = EmbeddedCurveParameters::AFFINE_GENERATOR_COEFFS;
        let generator: GroupAffine<EmbeddedCurveParameters> =
            GroupAffine::new(x, y);
        let point_f_pi: GroupAffine<EmbeddedCurveParameters> =
            AffineCurve::mul(
                &generator,
                EmbeddedCurveParameters::ScalarField::from(2u64).into_repr(),
            )
            .into_affine();

        // Prover POV
        let (proof, pi) = {
            let mut circuit: TestCircuit<
                P::ScalarField,
                EmbeddedCurveParameters,
            > = TestCircuit {
                a: P::ScalarField::from(20u64),
                b: P::ScalarField::from(5u64),
                c: P::ScalarField::from(25u64),
                d: P::ScalarField::from(100u64),
                e: EmbeddedCurveParameters::ScalarField::from(2u64),
                f: point_f_pi,
            };

            cfg_if::cfg_if! {
                if #[cfg(feature = "trace")] {
                    // Test trace
                    let mut prover: Prover<F, P, PC> = Prover::new(b"Test");
                    circuit.gadget(prover.mut_cs())?;
                    prover.cs.check_circuit_satisfied();
                }
            }

            circuit.gen_proof(&pp, pk, b"Test")?
        };

        let verifier_data = VerifierData::new(vk, pi);

        // Test serialisation for verifier_data
        let mut verifier_data_bytes = Vec::new();
        verifier_data.serialize(&mut verifier_data_bytes).unwrap();

        let deserialized_verifier_data: VerifierData<P> =
            VerifierData::deserialize(verifier_data_bytes.as_slice()).unwrap();

        assert!(deserialized_verifier_data == verifier_data);

        // Verifier POV

        // TODO: non-ideal hack for a first functional version.
        assert!(verify_proof::<P>(
            &pp,
            verifier_data.key,
            &proof,
            &verifier_data.pi,
            b"Test",
        )
        .is_ok());

        Ok(())
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_full_on_Bls12_381() -> Result<(), Error> {
        test_full::<
            Bls12_381_KZG,
            ark_bls12_381::Fr,
            ark_ed_on_bls12_381::EdwardsParameters,
        >()
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_full_on_Bls12_381_ipa() -> Result<(), Error> {
        test_full::<
            Bls12_381_IPA,
            ark_bls12_381::Fr,
            ark_ed_on_bls12_381::EdwardsParameters,
        >()
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_full_on_Bls12_377() -> Result<(), Error> {
        test_full::<
            Bls12_377_KZG,
            ark_bls12_377::Fr,
            ark_ed_on_bls12_377::EdwardsParameters,
        >()
    }
    #[test]
    #[allow(non_snake_case)]
    fn test_full_on_Bls12_377_ipa() -> Result<(), Error> {
        test_full::<
            Bls12_377_IPA,
            ark_bls12_377::Fr,
            ark_ed_on_bls12_377::EdwardsParameters,
        >()
    }
}
