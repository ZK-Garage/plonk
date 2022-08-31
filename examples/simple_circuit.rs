// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! PLONK Example

use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
use ark_ec::models::twisted_edwards_extended::GroupAffine;
use ark_ec::{AffineCurve, ProjectiveCurve, TEModelParameters};
use ark_ed_on_bls12_381::{
    EdwardsParameters as JubJubParameters, Fr as JubJubScalar,
};
use ark_ff::PrimeField;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly_commit::{sonic_pc::SonicKZG10, PolynomialCommitment};
use plonk::error::to_pc_error;
use plonk_core::circuit::{verify_proof, Circuit};
use plonk_core::constraint_system::StandardComposer;
use plonk_core::error::Error;
use plonk_core::prelude::*;
use rand_core::OsRng;

fn main() -> Result<(), Error> {
    // Implements a circuit that checks:
    // 1) a + b = c where C is a PI
    // 2) a <= 2^6
    // 3) b <= 2^4
    // 4) a * b = d where D is a PI
    // 5) JubJub::GENERATOR * e(JubJubScalar) = f where F is a PI
    #[derive(derivative::Derivative)]
    #[derivative(Debug(bound = ""), Default(bound = ""))]
    pub struct TestCircuit<F, P>
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
    {
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
            let zero = composer.zero_var();

            // Make first constraint a + b = c (as public input)
            composer.arithmetic_gate(|gate| {
                gate.witness(a, b, Some(zero))
                    .add(F::one(), F::one())
                    .pi(-self.c)
            });

            // Check that a and b are in range
            composer.range_gate(a, 6);
            composer.range_gate(b, 4);
            // Make second constraint a * b = d
            composer.arithmetic_gate(|gate| {
                gate.witness(a, b, Some(zero)).mul(F::one()).pi(-self.d)
            });
            let e =
                composer.add_input(from_embedded_curve_scalar::<F, P>(self.e));
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

    // Generate CRS
    type PC = SonicKZG10<Bls12_381, DensePolynomial<BlsScalar>>;
    let pp = PC::setup(1 << 10, None, &mut OsRng)
        .map_err(to_pc_error::<BlsScalar, PC>)?;

    let mut circuit = TestCircuit::<BlsScalar, JubJubParameters>::default();
    // Compile the circuit
    let (pk_p, (vk, _pi_pos)) = circuit.compile::<PC>(&pp)?;

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
        circuit.gen_proof::<PC>(&pp, pk_p, b"Test")
    }?;

    // Verifier POV
    let verifier_data = VerifierData::new(vk, pi);
    verify_proof::<BlsScalar, JubJubParameters, PC>(
        &pp,
        verifier_data.key,
        &proof,
        &verifier_data.pi,
        b"Test",
    )
}
