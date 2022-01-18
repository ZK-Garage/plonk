// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Tools & traits for PLONK circuits

use crate::{
    commitment::HomomorphicCommitment,
    error::{to_pc_error, Error},
    prelude::StandardComposer,
    proof_system::{Proof, Prover, ProverKey, Verifier, VerifierKey},
};
use ark_ec::models::TEModelParameters;
use ark_ff::{Field, PrimeField, ToConstraintField};
use ark_serialize::*;

/// Public Input Builder
#[derive(derivative::Derivative)]
#[derivative(Default(bound = ""))]
pub struct PublicInputBuilder<F>
where
    F: PrimeField,
{
    input: Vec<F>,
}

impl<F> PublicInputBuilder<F>
where
    F: PrimeField,
{
    /// Creates a new PublicInputBuilder
    pub fn new() -> Self {
        Self { input: vec![] }
    }
    /// Pushes a new input `item` into `self`, returning `None` if `item` could
    /// not be converted into constraint field elements.
    pub fn add_input<T>(self, item: &T) -> Option<Self>
    where
        T: ToConstraintField<F>,
    {
        self.extend([item])
    }
    /// Extends the input buffer by appending all the elements of `iter`
    /// converted into constraint field elements, returning `None` if any of
    /// the conversions failed.
    pub fn extend<'t, T, I>(mut self, iter: I) -> Option<Self>
    where
        T: 't + ToConstraintField<F>,
        I: IntoIterator<Item = &'t T>,
    {
        for item in iter {
            self.input.append(&mut item.to_field_elements()?);
        }
        Some(self)
    }
    /// Returns the vector of input elements.
    pub fn finish(self) -> Vec<F> {
        self.input
    }
}

/// Collection of structs/objects that the Verifier will use in order to
/// de/serialize data needed for Circuit proof verification.
/// This structure can be seen as a link between the [`Circuit`] public input
/// positions and the [`VerifierKey`] that the Verifier needs to use.
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = "VerifierKey<F,PC>: std::fmt::Debug"),
    Eq(bound = "VerifierKey<F,PC>: Eq"),
    PartialEq(bound = "VerifierKey<F,PC>: PartialEq")
)]
pub struct VerifierData<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    /// Verifier Key
    pub key: VerifierKey<F, PC>,

    /// Public Input Positions
    pub pi_pos: Vec<usize>,
}

impl<F, PC> VerifierData<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    /// Creates a new `VerifierData` from a [`VerifierKey`] and the public
    /// input positions of the circuit that it represents.
    pub fn new(key: VerifierKey<F, PC>, pi_pos: Vec<usize>) -> Self {
        Self { key, pi_pos }
    }

    /// Returns a reference to the contained [`VerifierKey`].
    pub fn key(&self) -> &VerifierKey<F, PC> {
        &self.key
    }

    /// Returns a reference to the contained Public Input positions.
    pub fn pi_pos(&self) -> &[usize] {
        &self.pi_pos
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
/// use plonk_core::circuit::{Circuit, verify_proof, PublicInputBuilder};
/// use plonk_core::constraint_system::StandardComposer;
/// use plonk_core::error::{to_pc_error,Error};
/// use ark_poly::polynomial::univariate::DensePolynomial;
/// use ark_poly_commit::{PolynomialCommitment, sonic_pc::SonicKZG10};
/// use plonk_core::prelude::*;
/// use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
/// use num_traits::{Zero, One};
/// use rand::rngs::OsRng;
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
/// pub struct TestCircuit<F, P>
/// where
///     F: PrimeField,
///     P: TEModelParameters<BaseField = F>,
/// {
///        a: F,
///        b: F,
///        c: F,
///        d: F,
///        e: P::ScalarField,
///        f: GroupAffine<P>,
///    }
///
/// impl<F, P> Circuit<F, P> for TestCircuit<F, P>
/// where
///     F: PrimeField,
///     P: TEModelParameters<BaseField = F>,
///    {
///        const CIRCUIT_ID: [u8; 32] = [0xff; 32];
///
///        fn gadget(
///            &mut self,
///            composer: &mut StandardComposer<F, P>,
///        ) -> Result<(), Error> {
///            let a = composer.add_input(self.a);
///            let b = composer.add_input(self.b);
///            let zero = composer.zero_var();
///
///            // Make first constraint a + b = c (as public input)
///            composer.arithmetic_gate(|gate| {
///                gate.witness(a, b, Some(zero))
///                    .add(F::one(), F::one())
///                    .pi(-self.c)
///            });
///
///            // Check that a and b are in range
///            composer.range_gate(a, 1 << 6);
///            composer.range_gate(b, 1 << 5);
///            // Make second constraint a * b = d
///            composer.arithmetic_gate(|gate| {
///                gate.witness(a, b, Some(zero)).mul(F::one()).pi(-self.d)
///            });
///            let e = composer
///                .add_input(from_embedded_curve_scalar::<F, P>(self.e));
///            let (x, y) = P::AFFINE_GENERATOR_COEFFS;
///            let generator = GroupAffine::new(x, y);
///            let scalar_mul_result =
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
/// let pp = PC::setup(
///     1 << 12, None, &mut OsRng
///  )?;
///
/// let mut circuit = TestCircuit::<BlsScalar, JubJubParameters>::default();
/// // Compile the circuit
/// let (pk_p, verifier_data) = circuit.compile::<PC>(&pp)?;
///
/// let (x, y) = JubJubParameters::AFFINE_GENERATOR_COEFFS;
/// let generator: GroupAffine<JubJubParameters> = GroupAffine::new(x, y);
/// let point_f_pi: GroupAffine<JubJubParameters> = AffineCurve::mul(
///     &generator,
///     JubJubScalar::from(2u64).into_repr(),
/// )
/// .into_affine();
/// // Prover POV
/// let proof = {
///     let mut circuit: TestCircuit<BlsScalar, JubJubParameters> = TestCircuit {
///         a: BlsScalar::from(20u64),
///         b: BlsScalar::from(5u64),
///         c: BlsScalar::from(25u64),
///         d: BlsScalar::from(100u64),
///         e: JubJubScalar::from(2u64),
///         f: point_f_pi,
///     };
///     circuit.gen_proof::<PC>(&pp, pk_p, b"Test")
/// }?;
///
/// // Test serialisation for verifier_data
/// let mut verifier_data_bytes = Vec::new();
/// verifier_data.serialize(&mut verifier_data_bytes).unwrap();
///
/// let deserialized_verifier_data: VerifierData<BlsScalar, PC> =
///     VerifierData::deserialize(verifier_data_bytes.as_slice()).unwrap();
///
/// assert!(deserialized_verifier_data == verifier_data);
///
/// // Verifier POV
/// let public_inputs = PublicInputBuilder::new()
///     .add_input(&BlsScalar::from(25u64))
///     .unwrap()
///     .add_input(&BlsScalar::from(100u64))
///     .unwrap()
///     .add_input(&point_f_pi)
///     .unwrap()
///     .finish();
///
/// let VerifierData { key, pi_pos } = verifier_data;
/// // TODO: non-ideal hack for a first functional version.
/// verify_proof::<BlsScalar, JubJubParameters, PC>(
///     &pp,
///     key,
///     &proof,
///     &public_inputs,
///     &pi_pos,
///     b"Test",
/// )
/// }
/// ```
pub trait Circuit<F, P>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    /// Circuit identifier associated constant.
    const CIRCUIT_ID: [u8; 32];

    /// Gadget implementation used to fill the composer.
    fn gadget(
        &mut self,
        composer: &mut StandardComposer<F, P>,
    ) -> Result<(), Error>;

    /// Compiles the circuit by using a function that returns a `Result`
    /// with the `ProverKey`, `VerifierKey` and the circuit size.
    #[allow(clippy::type_complexity)] // NOTE: Clippy is too hash here.
    fn compile<PC>(
        &mut self,
        u_params: &PC::UniversalParams,
    ) -> Result<(ProverKey<F>, VerifierData<F, PC>), Error>
    where
        F: PrimeField,
        PC: HomomorphicCommitment<F>,
    {
        // Setup PublicParams
        let circuit_size = self.padded_circuit_size();
        let (ck, _) = PC::trim(
            u_params,
            // +1 per wire, +2 for the permutation poly
            circuit_size + 6,
            0,
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;

        //Generate & save `ProverKey` with some random values.
        let mut prover = Prover::<F, P, PC>::new(b"CircuitCompilation");
        self.gadget(prover.mut_cs())?;
        let pi_pos = prover.mut_cs().pi_positions();
        prover.preprocess(&ck)?;

        // Generate & save `VerifierKey` with some random values.
        let mut verifier = Verifier::new(b"CircuitCompilation");
        self.gadget(verifier.mut_cs())?;
        verifier.preprocess(&ck)?;
        Ok((
            prover
                .prover_key
                .expect("Unexpected error. Missing ProverKey in compilation"),
            VerifierData::new(
                verifier.verifier_key.expect(
                    "Unexpected error. Missing VerifierKey in compilation",
                ),
                pi_pos,
            ),
        ))
    }

    /// Generates a proof using the provided `CircuitInputs` & `ProverKey`
    /// instances.
    fn gen_proof<PC>(
        &mut self,
        u_params: &PC::UniversalParams,
        prover_key: ProverKey<F>,
        transcript_init: &'static [u8],
    ) -> Result<Proof<F, PC>, Error>
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let circuit_size = self.padded_circuit_size();
        let (ck, _) = PC::trim(
            u_params,
            // +1 per wire, +2 for the permutation poly
            circuit_size + 6,
            0,
            None,
        )
        .map_err(to_pc_error::<F, PC>)?;
        // New Prover instance
        let mut prover = Prover::new(transcript_init);
        // Fill witnesses for Prover
        self.gadget(prover.mut_cs())?;
        // Add ProverKey to Prover
        prover.prover_key = Some(prover_key);
        prover.prove(&ck)
    }

    /// Returns the Circuit size padded to the next power of two.
    fn padded_circuit_size(&self) -> usize;
}

/// Verifies a proof using the provided `CircuitInputs` & `VerifierKey`
/// instances.
pub fn verify_proof<F, P, PC>(
    u_params: &PC::UniversalParams,
    plonk_verifier_key: VerifierKey<F, PC>,
    proof: &Proof<F, PC>,
    pub_inputs_values: &[F],
    pub_inputs_positions: &[usize],
    transcript_init: &'static [u8],
) -> Result<(), Error>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
    PC: HomomorphicCommitment<F>,
{
    let mut verifier: Verifier<F, P, PC> = Verifier::new(transcript_init);
    let padded_circuit_size = plonk_verifier_key.padded_circuit_size();
    verifier.verifier_key = Some(plonk_verifier_key);
    let (_, vk) = PC::trim(
        u_params,
        // +1 per wire, +2 for the permutation poly
        padded_circuit_size + 6,
        0,
        None,
    )
    .map_err(to_pc_error::<F, PC>)?;

    verifier.verify(
        proof,
        &vk,
        build_pi(pub_inputs_values, pub_inputs_positions, padded_circuit_size)
            .as_slice(),
    )
}

/// Build PI vector for Proof verifications.
fn build_pi<'a, F>(
    pub_input_values: impl IntoIterator<Item = &'a F>,
    pub_input_pos: &[usize],
    trim_size: usize,
) -> Vec<F>
where
    F: Field,
{
    let mut pi = vec![F::zero(); trim_size];
    pub_input_values
        .into_iter()
        .zip(pub_input_pos.iter().copied())
        .for_each(|(value, pos)| {
            pi[pos] = -*value;
        });
    pi
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{constraint_system::StandardComposer, util};
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_ec::{
        twisted_edwards_extended::GroupAffine, AffineCurve, PairingEngine,
        ProjectiveCurve,
    };
    use ark_ff::{FftField, PrimeField};
    use rand::rngs::OsRng;

    // Implements a circuit that checks:
    // 1) a + b = c where C is a PI
    // 2) a <= 2^6
    // 3) b <= 2^5
    // 4) a * b = d where D is a PI
    // 5) JubJub::GENERATOR * e(JubJubScalar) = f where F is a PI
    #[derive(derivative::Derivative)]
    #[derivative(Debug(bound = ""), Default(bound = ""))]
    pub struct TestCircuit<F: FftField, P: TEModelParameters<BaseField = F>> {
        a: F,
        b: F,
        c: F,
        d: F,
        e: P::ScalarField,
        f: GroupAffine<P>,
    }

    impl<F, P> Circuit<F, P> for TestCircuit<F, P>
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
    {
        const CIRCUIT_ID: [u8; 32] = [0xff; 32];

        fn gadget(
            &mut self,
            composer: &mut StandardComposer<F, P>,
        ) -> Result<(), Error> {
            let a = composer.add_input(self.a);
            let b = composer.add_input(self.b);
            let zero = composer.zero_var;

            // Make first constraint a + b = c (as public input)
            composer.arithmetic_gate(|gate| {
                gate.witness(a, b, Some(zero))
                    .add(F::one(), F::one())
                    .pi(-self.c)
            });

            // Check that a and b are in range
            composer.range_gate(a, 1 << 6);
            composer.range_gate(b, 1 << 5);
            // Make second constraint a * b = d
            composer.arithmetic_gate(|gate| {
                gate.witness(a, b, Some(zero)).mul(F::one()).pi(-self.d)
            });
            let e = composer
                .add_input(util::from_embedded_curve_scalar::<F, P>(self.e));
            let (x, y) = P::AFFINE_GENERATOR_COEFFS;
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

    fn test_full<F, P, PC>() -> Result<(), Error>
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
        VerifierData<F, PC>: PartialEq,
    {
        // Generate CRS
        let pp = PC::setup(1 << 19, None, &mut OsRng)
            .map_err(to_pc_error::<F, PC>)?;

        let mut circuit = TestCircuit::<F, P>::default();

        // Compile the circuit
        let (pk_p, verifier_data) = circuit.compile::<PC>(&pp)?;

        let (x, y) = P::AFFINE_GENERATOR_COEFFS;
        let generator: GroupAffine<P> = GroupAffine::new(x, y);
        let point_f_pi: GroupAffine<P> = AffineCurve::mul(
            &generator,
            P::ScalarField::from(2u64).into_repr(),
        )
        .into_affine();

        // Prover POV
        let proof = {
            let mut circuit: TestCircuit<F, P> = TestCircuit {
                a: F::from(20u64),
                b: F::from(5u64),
                c: F::from(25u64),
                d: F::from(100u64),
                e: P::ScalarField::from(2u64),
                f: point_f_pi,
            };

            circuit.gen_proof::<PC>(&pp, pk_p, b"Test")?
        };

        // Test serialisation for verifier_data
        let mut verifier_data_bytes = Vec::new();
        verifier_data.serialize(&mut verifier_data_bytes).unwrap();

        let deserialized_verifier_data: VerifierData<F, PC> =
            VerifierData::deserialize(verifier_data_bytes.as_slice()).unwrap();

        assert!(deserialized_verifier_data == verifier_data);

        // Verifier POV
        let public_inputs = PublicInputBuilder::new()
            .add_input(&F::from(25u64))
            .unwrap()
            .add_input(&F::from(100u64))
            .unwrap()
            .add_input(&point_f_pi)
            .unwrap()
            .finish();

        let VerifierData { key, pi_pos } = verifier_data;

        // TODO: non-ideal hack for a first functional version.
        assert!(verify_proof::<F, P, PC>(
            &pp,
            key,
            &proof,
            &public_inputs,
            &pi_pos,
            b"Test",
        )
        .is_ok());

        Ok(())
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_full_on_Bls12_381() -> Result<(), Error> {
        test_full::<
            <Bls12_381 as PairingEngine>::Fr,
            ark_ed_on_bls12_381::EdwardsParameters,
            crate::commitment::KZG10<Bls12_381>,
        >()
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_full_on_Bls12_381_ipa() -> Result<(), Error> {
        test_full::<
            <Bls12_381 as PairingEngine>::Fr,
            ark_ed_on_bls12_381::EdwardsParameters,
            crate::commitment::IPA<
                <Bls12_381 as PairingEngine>::G1Affine,
                blake2::Blake2b,
            >,
        >()
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_full_on_Bls12_377() -> Result<(), Error> {
        test_full::<
            <Bls12_377 as PairingEngine>::Fr,
            ark_ed_on_bls12_377::EdwardsParameters,
            crate::commitment::KZG10<Bls12_377>,
        >()
    }
    #[test]
    #[allow(non_snake_case)]
    fn test_full_on_Bls12_377_ipa() -> Result<(), Error> {
        test_full::<
            <Bls12_377 as PairingEngine>::Fr,
            ark_ed_on_bls12_377::EdwardsParameters,
            crate::commitment::IPA<
                <Bls12_377 as PairingEngine>::G1Affine,
                blake2::Blake2b,
            >,
        >()
    }
}
