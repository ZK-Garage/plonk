// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Tools & traits for PLONK circuits

use crate::constraint_system::StandardComposer;
use crate::error::Error;
use crate::proof_system::{Proof, Prover, ProverKey, Verifier, VerifierKey};
use ark_ec::models::{SWModelParameters, TEModelParameters};
use ark_ec::{
    short_weierstrass_jacobian::{
        GroupAffine as SWGroupAffine, GroupProjective as SWGroupProjective,
    },
    twisted_edwards_extended::{
        GroupAffine as TEGroupAffine, GroupProjective as TEGroupProjective,
    },
    PairingEngine, ProjectiveCurve,
};
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::kzg10::{self, Powers, UniversalParams};
use ark_poly_commit::sonic_pc::SonicKZG10;
use ark_poly_commit::PolynomialCommitment;
use ark_serialize::*;

/// Field Element Into Public Input
///
/// The reason for introducing these two traits, `FeIntoPubInput` and
/// `GeIntoPubInput` is to have a workaround for not being able to
/// implement `From<_> for Values` for both `PrimeField` and `GroupAffine`. The
/// reason why this is not possible is because both the trait `PrimeField` and
/// the struct `GroupAffine` are external to the crate, and therefore the
/// compiler cannot be sure that `PrimeField` will never be implemented for
/// `GroupAffine`. In which case, the two implementations of `From` would be
/// inconsistent. To this end, we create to helper traits, `FeIntoPubInput` and
/// `GeIntoPubInput`, that stand for "Field Element Into Public Input" and
/// "Group Element Into Public Input" respectively.
pub trait FeIntoPubInput<T> {
    /// Ad hoc `Into` implementation. Serves the same purpose as `Into`, but as
    /// a different trait. Read documentation of Trait for more details.
    fn into_pi(self) -> T;
}

/// Group Element Into Public Input
///
/// The reason for introducing these two traits is to have a workaround for not
/// being able to implement `From<_> for Values` for both `PrimeField` and
/// `GroupAffine`. The reason why this is not possible is because both the trait
/// `PrimeField` and the struct `GroupAffine` are external to the crate, and
/// therefore the compiler cannot be sure that `PrimeField` will never be
/// implemented for `GroupAffine`. In which case, the two implementations of
/// `From` would be inconsistent. To this end, we create to helper traits,
/// `FeIntoPubInput` and `GeIntoPubInput`, that stand for "Field Element Into
/// Public Input" and "Group Element Into Public Input" respectively.
pub trait GeIntoPubInput<T> {
    /// Ad hoc `Into` implementation. Serves the same purpose as `Into`, but as
    /// a different trait. Read documentation of Trait for more details.
    fn into_pi(self) -> T;
}

/// Structure that represents a PLONK Circuit Public Input converted into its
/// scalar representation.
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Default)]
pub struct PublicInputValue<F>
where
    F: Field,
{
    /// Field Values
    pub(crate) values: Vec<F>,
}
/*
impl<F> FeIntoPubInput<PublicInputValue<F>> for F
where
    F : Field,
{
    #[inline]
    fn into_pi(self) -> PublicInputValue<F> {
        PublicInputValue { values: vec![self] }
    }
}*/

impl<F> From<F> for PublicInputValue<F>
where
    F: Field,
{
    fn from(p: F) -> PublicInputValue<F> {
        PublicInputValue { values: vec![p] }
    }
}

impl<P> GeIntoPubInput<PublicInputValue<P::BaseField>> for TEGroupAffine<P>
where
    P: TEModelParameters,
{
    #[inline]
    fn into_pi(self) -> PublicInputValue<P::BaseField> {
        PublicInputValue {
            values: vec![self.x, self.y],
        }
    }
}

impl<P> GeIntoPubInput<PublicInputValue<P::BaseField>> for TEGroupProjective<P>
where
    P: TEModelParameters,
{
    #[inline]
    fn into_pi(self) -> PublicInputValue<P::BaseField> {
        GeIntoPubInput::into_pi(self.into_affine())
    }
}
impl<P> GeIntoPubInput<PublicInputValue<P::BaseField>> for SWGroupAffine<P>
where
    P: SWModelParameters,
{
    #[inline]
    fn into_pi(self) -> PublicInputValue<P::BaseField> {
        PublicInputValue {
            values: vec![self.x, self.y],
        }
    }
}

impl<P> GeIntoPubInput<PublicInputValue<P::BaseField>> for SWGroupProjective<P>
where
    P: SWModelParameters,
{
    #[inline]
    fn into_pi(self) -> PublicInputValue<P::BaseField> {
        GeIntoPubInput::into_pi(self.into_affine())
    }
}

/// Collection of structs/objects that the Verifier will use in order to
/// de/serialize data needed for Circuit proof verification.
/// This structure can be seen as a link between the [`Circuit`] public input
/// positions and the [`VerifierKey`] that the Verifier needs to use.
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
    Eq(bound = ""),
    PartialEq(bound = "")
)]
pub struct VerifierData<E>
where
    E: PairingEngine,
{
    /// Verifier Key
    pub key: VerifierKey<E>,

    /// Public Input Positions
    pub pi_pos: Vec<usize>,
}

impl<E> VerifierData<E>
where
    E: PairingEngine,
{
    /// Creates a new `VerifierData` from a [`VerifierKey`] and the public
    /// input positions of the circuit that it represents.
    pub fn new(key: VerifierKey<E>, pi_pos: Vec<usize>) -> Self {
        Self { key, pi_pos }
    }

    /// Returns a reference to the contained [`VerifierKey`].
    pub fn key(&self) -> &VerifierKey<E> {
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
/// ```rust
/// use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
/// use ark_ec::PairingEngine;
/// use ark_ec::models::twisted_edwards_extended::GroupAffine;
/// use ark_ec::{TEModelParameters, AffineCurve, ProjectiveCurve};
/// use ark_ed_on_bls12_381::{
///     EdwardsAffine as JubjubAffine, EdwardsParameters as JubjubParameters,
///     EdwardsProjective as JubjubProjective, Fr as JubjubScalar,
/// };
/// use ark_ff::{FftField, PrimeField, BigInteger};
/// use ark_plonk::circuit::{Circuit, PublicInputValue, verify_proof, GeIntoPubInput, FeIntoPubInput};
/// use ark_plonk::constraint_system::StandardComposer;
/// use ark_plonk::error::Error;
/// use ark_plonk::prelude::VerifierData;
/// use ark_poly::polynomial::univariate::DensePolynomial;
/// use ark_poly_commit::kzg10::KZG10;
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
/// #[derive(Debug)]
/// pub struct TestCircuit<F, P>
/// where
///     F: FftField,
///     P: TEModelParameters<BaseField = F>,
/// {
///     a: F,
///     b: F,
///     c: F,
///     d: F,
///     e: P::ScalarField,
///     f: GroupAffine<P>,
/// }
///
/// impl<F, P> Default for TestCircuit<F, P>
/// where
///     F: FftField,
///     P: TEModelParameters<BaseField = F>
/// {
///     fn default() -> Self {
///         Self {
///             a: F::zero(),
///             b: F::zero(),
///             c: F::zero(),
///             d: F::zero(),
///             e: P::ScalarField::zero(),
///             f: GroupAffine::<P>::zero(),
///         }
///     }
/// }
///
/// impl<F, P> Circuit<F> for TestCircuit<F, P>
/// where
///     F: FftField + PrimeField,
///     P: TEModelParameters<BaseField = F>,
/// {
///     const CIRCUIT_ID: [u8; 32] = [0xff; 32];
///
///     fn gadget(
///         &mut self,
///         composer: &mut StandardComposer<F>,
///     ) -> Result<(), Error> {
///         // Add fixed witness zero
///         let a = composer.add_input(self.a);
///         let b = composer.add_input(self.b);
///         // Make first constraint a + b = c
///         let add_result = composer.add(
///           (F::one(), a),
///           (F::one(), b),
///           F::zero(),
///           Some(-self.c),
///         );
///         composer.assert_equal(add_result, composer.zero_var());
///
///         // Check that a and b are in range
///         composer.range_gate(a, 1 << 6);
///         composer.range_gate(b, 1 << 5);
///
///         // Make second constraint a * b = d
///         let mul_result = composer.mul(F::one(), a, b, F::zero(), Some(-self.d));
///         composer.assert_equal(mul_result, composer.zero_var());
///
///         let e_repr = self.e.into_repr().to_bytes_le();
///         let e = composer.add_input(F::from_le_bytes_mod_order(&e_repr));
///         let (x, y) = P::AFFINE_GENERATOR_COEFFS;
///         let generator = GroupAffine::<P>::new(x, y);
///         let scalar_mul_result =
///             composer.fixed_base_scalar_mul(e, generator);
///         // Apply the constraint
///         composer
///             .assert_equal_public_point(scalar_mul_result, self.f.clone());
///         Ok(())
///     }
///
///     fn padded_circuit_size(&self) -> usize {
///         1 << 11
///     }
/// }
///
/// let pp = KZG10::<Bls12_381,DensePolynomial<BlsScalar>,>::setup(
///     1 << 12, false, &mut OsRng
///  )?;
///
/// // Initialize the circuit
/// let mut circuit = TestCircuit::<BlsScalar, JubjubParameters>::default();
///
/// // Compile the circuit
/// let (pk, vd) = circuit.compile::<Bls12_381>(&pp)?;
///
/// // Prover POV
/// let (x, y) = JubjubParameters::AFFINE_GENERATOR_COEFFS;
/// let generator = JubjubAffine::new(x, y);
/// let point_f_pi: JubjubAffine = AffineCurve::mul(
///   &generator,
///   JubjubScalar::from(2u64).into_repr(),
/// ).into_affine();
///
/// let proof = {
///     let mut circuit = TestCircuit {
///         a: BlsScalar::from(20u64),
///         b: BlsScalar::from(5u64),
///         c: BlsScalar::from(25u64),
///         d: BlsScalar::from(100u64),
///         e: JubjubScalar::from(2u64),
///         f: point_f_pi,
///     };
///     circuit.gen_proof::<Bls12_381, JubjubParameters>(&pp, pk, b"Test")
/// }?;
///
/// // Verifier POV
/// let public_inputs: Vec<PublicInputValue<BlsScalar>> = vec![
///     BlsScalar::from(25u64).into(),
///     BlsScalar::from(100u64).into(),
///     GeIntoPubInput::into_pi(point_f_pi),
/// ];
/// let VerifierData { key, pi_pos } = vd;
/// verify_proof::<Bls12_381,JubjubParameters>(
///     &pp,
///     key,
///     &proof,
///     &public_inputs,
///     &pi_pos,
///     b"Test",
/// )
/// }
/// ```
pub trait Circuit<F>
where
    F: FftField,
{
    /// Circuit identifier associated constant.
    const CIRCUIT_ID: [u8; 32];

    /// Gadget implementation used to fill the composer.
    fn gadget(
        &mut self,
        composer: &mut StandardComposer<F>,
    ) -> Result<(), Error>;

    /// Compiles the circuit by using a function that returns a `Result`
    /// with the `ProverKey`, `VerifierKey` and the circuit size.
    #[allow(clippy::type_complexity)] // NOTE: Clippy is too hash here.
    fn compile<E: PairingEngine<Fr = F>>(
        &mut self,
        u_params: &UniversalParams<E>,
    ) -> Result<(ProverKey<E::Fr>, VerifierData<E>), Error>
    where
        F: PrimeField,
    {
        // Setup PublicParams
        // XXX: KZG10 does not have a trim function so we use sonics and
        // then do a transformation between sonic CommiterKey to KZG10
        // powers
        let circuit_size = self.padded_circuit_size();
        let (ck, _) = SonicKZG10::<E, DensePolynomial<E::Fr>>::trim(
            u_params,
            circuit_size,
            0,
            None,
        )
        .unwrap();
        let powers = Powers {
            powers_of_g: ck.powers_of_g.into(),
            powers_of_gamma_g: ck.powers_of_gamma_g.into(),
        };
        //Generate & save `ProverKey` with some random values.
        let mut prover = Prover::new(b"CircuitCompilation");
        self.gadget(prover.mut_cs())?;
        let pi_pos = prover.mut_cs().pi_positions();
        prover.preprocess(&powers)?;

        // Generate & save `VerifierKey` with some random values.
        let mut verifier = Verifier::new(b"CircuitCompilation");
        self.gadget(verifier.mut_cs())?;
        verifier.preprocess(&powers)?;
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
    fn gen_proof<E, P>(
        &mut self,
        u_params: &UniversalParams<E>,
        prover_key: ProverKey<F>,
        transcript_init: &'static [u8],
    ) -> Result<Proof<E>, Error>
    where
        F: PrimeField,
        E: PairingEngine<Fr = F>,
        P: TEModelParameters<BaseField = F>,
    {
        // XXX: KZG10 does not have a trim function so we use sonics and
        // then do a transformation between sonic CommiterKey to KZG10
        // powers
        let circuit_size = self.padded_circuit_size();
        let (ck, _) = SonicKZG10::<E, DensePolynomial<E::Fr>>::trim(
            u_params,
            circuit_size,
            0,
            None,
        )
        .unwrap();
        let powers = Powers::<E> {
            powers_of_g: ck.powers_of_g.into(),
            powers_of_gamma_g: ck.powers_of_gamma_g.into(),
        };
        // New Prover instance
        let mut prover = Prover::new(transcript_init);
        // Fill witnesses for Prover
        self.gadget(prover.mut_cs())?;
        // Add ProverKey to Prover
        prover.prover_key = Some(prover_key);
        prover.prove::<E, P>(&powers)
    }

    /// Returns the Circuit size padded to the next power of two.
    fn padded_circuit_size(&self) -> usize;
}

/// Verifies a proof using the provided `CircuitInputs` & `VerifierKey`
/// instances.
pub fn verify_proof<E, P>(
    u_params: &UniversalParams<E>,
    plonk_verifier_key: VerifierKey<E>,
    proof: &Proof<E>,
    pub_inputs_values: &[PublicInputValue<E::Fr>],
    pub_inputs_positions: &[usize],
    transcript_init: &'static [u8],
) -> Result<(), Error>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    let mut verifier: Verifier<E> = Verifier::new(transcript_init);
    let padded_circuit_size = plonk_verifier_key.padded_circuit_size();
    // let key: VerifierKey<E, P> = *plonk_verifier_key;
    verifier.verifier_key = Some(plonk_verifier_key);
    let (_, sonic_vk) = SonicKZG10::<E, DensePolynomial<E::Fr>>::trim(
        u_params,
        padded_circuit_size,
        0,
        None,
    )
    .unwrap();

    let vk = kzg10::VerifierKey {
        g: sonic_vk.g,
        gamma_g: sonic_vk.gamma_g,
        h: sonic_vk.h,
        beta_h: sonic_vk.beta_h,
        prepared_h: sonic_vk.prepared_h,
        prepared_beta_h: sonic_vk.prepared_beta_h,
    };

    verifier.verify::<P>(
        proof,
        &vk,
        build_pi(pub_inputs_values, pub_inputs_positions, padded_circuit_size)
            .as_slice(),
    )
}

/// Build PI vector for Proof verifications.
fn build_pi<F>(
    pub_input_values: &[PublicInputValue<F>],
    pub_input_pos: &[usize],
    trim_size: usize,
) -> Vec<F>
where
    F: Field,
{
    let mut pi = vec![F::zero(); trim_size];
    pub_input_values
        .iter()
        .flat_map(|pub_input| pub_input.values.clone())
        .zip(pub_input_pos.iter().copied())
        .for_each(|(value, pos)| {
            pi[pos] = -value;
        });
    pi
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{constraint_system::StandardComposer, util};
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_ec::twisted_edwards_extended::GroupAffine;
    use ark_ec::AffineCurve;
    use ark_ff::{FftField, PrimeField};
    use ark_poly_commit::kzg10::KZG10;

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

    impl<F, P> Circuit<F> for TestCircuit<F, P>
    where
        F: FftField + PrimeField,
        P: TEModelParameters<BaseField = F>,
    {
        const CIRCUIT_ID: [u8; 32] = [0xff; 32];

        fn gadget(
            &mut self,
            composer: &mut StandardComposer<F>,
        ) -> Result<(), Error> {
            let a = composer.add_input(self.a);
            let b = composer.add_input(self.b);

            // Make first constraint a + b = c (as public input)
            let add_result = composer.add(
                (F::one(), a),
                (F::one(), b),
                F::zero(),
                Some(-self.c),
            );
            composer.assert_equal(add_result, composer.zero_var());

            // Check that a and b are in range
            composer.range_gate(a, 1 << 6);
            composer.range_gate(b, 1 << 5);
            // Make second constraint a * b = d
            let mul_result =
                composer.mul(F::one(), a, b, F::zero(), Some(-self.d));
            composer.assert_equal(mul_result, composer.zero_var());

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
            1 << 11
        }
    }

    fn test_full<E: PairingEngine, P: TEModelParameters<BaseField = E::Fr>>(
    ) -> Result<(), Error> {
        use rand_core::OsRng;

        // Generate CRS
        let pp = KZG10::<E, DensePolynomial<E::Fr>>::setup(
            1 << 12,
            false,
            &mut OsRng,
        )?;

        let mut circuit = TestCircuit::<E::Fr, P>::default();

        // Compile the circuit
        let (pk_p, verifier_data) = circuit.compile(&pp)?;

        let (x, y) = P::AFFINE_GENERATOR_COEFFS;
        let generator: GroupAffine<P> = GroupAffine::new(x, y);
        let point_f_pi: GroupAffine<P> = AffineCurve::mul(
            &generator,
            P::ScalarField::from(2u64).into_repr(),
        )
        .into_affine();

        // Prover POV
        let proof = {
            let mut circuit: TestCircuit<E::Fr, P> = TestCircuit {
                a: E::Fr::from(20u64),
                b: E::Fr::from(5u64),
                c: E::Fr::from(25u64),
                d: E::Fr::from(100u64),
                e: P::ScalarField::from(2u64),
                f: point_f_pi,
            };

            circuit.gen_proof::<E, P>(&pp, pk_p, b"Test")?
        };

        // Test serialisation for verifier_data
        let mut verifier_data_bytes = Vec::new();
        verifier_data.serialize(&mut verifier_data_bytes).unwrap();

        let verif_data: VerifierData<E> =
            VerifierData::deserialize(verifier_data_bytes.as_slice()).unwrap();

        assert!(verif_data == verifier_data);

        // Verifier POV
        let public_inputs: Vec<PublicInputValue<E::Fr>> = vec![
            E::Fr::from(25u64).into(),
            E::Fr::from(100u64).into(),
            GeIntoPubInput::into_pi(point_f_pi),
        ];

        let VerifierData { key, pi_pos } = verifier_data;

        // TODO: non-ideal hack for a first functional version.
        assert!(verify_proof::<E, P>(
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
        test_full::<Bls12_381, ark_ed_on_bls12_381::EdwardsParameters>()
    }

    #[test]
    #[allow(non_snake_case)]
    fn test_full_on_Bls12_377() -> Result<(), Error> {
        test_full::<Bls12_377, ark_ed_on_bls12_377::EdwardsParameters>()
    }
}
