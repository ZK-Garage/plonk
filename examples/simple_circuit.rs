// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! PLONK Example

use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
use ark_ec::{
    models::twisted_edwards_extended::GroupAffine, AffineCurve,
    ProjectiveCurve, TEModelParameters,
};
use ark_ed_on_bls12_381::{
    EdwardsParameters as JubJubParameters, Fr as JubJubScalar,
};
use ark_ff::{FftField, PrimeField};
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly_commit::{sonic_pc::SonicKZG10, PolynomialCommitment};
use plonk::error::to_pc_error;
use plonk_core::{
    circuit::{verify_proof, Circuit},
    constraint_system::StandardComposer,
    error::Error,
    prelude::*,
};
use rand_core::OsRng;

fn main() -> Result<(), Error> {
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
            let zero = composer.zero_var();

            // Make first constraint a + b = c (as public input)
            composer.arithmetic_gate(|gate| {
                gate.witness(a, b, Some(zero))
                    .add(EmbeddedBaseField::one(), EmbeddedBaseField::one())
                    .pi(-self.c)
            });

            // Check that a and b are in range
            composer.range_gate(a, 1 << 6);
            composer.range_gate(b, 1 << 5);
            // Make second constraint a * b = d
            composer.arithmetic_gate(|gate| {
                gate.witness(a, b, Some(zero))
                    .mul(EmbeddedBaseField::one())
                    .pi(-self.d)
            });
            let e = composer.add_input(from_embedded_curve_scalar::<
                EmbeddedBaseField,
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

    // Generate CRS
    type PC = SonicKZG10<Bls12_381, DensePolynomial<BlsScalar>>;
    type P = plonk_core::parameters::test::Bls12_381_KZG;
    let pp = PC::setup(1 << 10, None, &mut OsRng).map_err(to_pc_error::<P>)?;

    let mut circuit = TestCircuit::<BlsScalar, JubJubParameters>::default();
    // Compile the circuit
    let (pk_p, vk) = circuit.compile(&pp)?;

    let (x, y) = JubJubParameters::AFFINE_GENERATOR_COEFFS;
    let generator: GroupAffine<JubJubParameters> = GroupAffine::new(x, y);
    let point_f_pi: GroupAffine<JubJubParameters> =
        AffineCurve::mul(&generator, JubJubScalar::from(2u64).into_repr())
            .into_affine();
    // Prover POV
    let (proof, pi) = {
        let mut circuit: TestCircuit<BlsScalar, JubJubParameters> =
            TestCircuit {
                a: BlsScalar::from(20u64),
                b: BlsScalar::from(5u64),
                c: BlsScalar::from(25u64),
                d: BlsScalar::from(100u64),
                e: JubJubScalar::from(2u64),
                f: point_f_pi,
            };
        circuit.gen_proof(&pp, pk_p, b"Test")
    }?;

    // Verifier POV
    let verifier_data = VerifierData::new(vk, pi);
    verify_proof::<P>(
        &pp,
        verifier_data.key,
        &proof,
        &verifier_data.pi,
        b"Test",
    )
}
