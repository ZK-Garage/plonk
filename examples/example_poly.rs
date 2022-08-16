// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! PLONK Example

use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
use ark_ec::TEModelParameters;
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
    // Implements a circuit that defines the following function:
    // -------------------
    // def f(x, y, z):
    //   if x == 1:
    //     return y*z
    //   return 2y - z
    // -------------------
    // where x,y,z are inputs and r is the output

    #[derive(derivative::Derivative)]
    #[derivative(Debug(bound = ""), Default(bound = ""))]
    pub struct TestCircuit<F, P>
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
    {
        x: F,
        y: F,
        z: F,
        r: F,
        dummy: P::ScalarField,
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
            let x = composer.add_input(self.x);
            let y = composer.add_input(self.y);
            let z = composer.add_input(self.z);
            let r = composer.add_input(self.r);

            let one = composer.add_input(F::one());
            let zero = composer.zero_var();

            // Define 2y-z
            let two_y_minus_z = composer.arithmetic_gate(|gate| {
                gate.witness(y, z, None).add(F::from(2u64), -F::one())
            });

            // Define y*z
            let y_times_z = composer
                .arithmetic_gate(|gate| gate.witness(y, z, None).mul(F::one()));

            // Define x-1
            let x_minus_1 = composer.arithmetic_gate(|gate| {
                gate.witness(x, zero, None)
                    .add(F::one(), F::zero())
                    .constant(-F::one())
            });

            // Define x-1==0
            let x_bool = composer.is_zero_with_output(x_minus_1);

            // Define x bool negate
            // x_bool_negate also be implemented with an add gate as
            // x_bool_negate = 1 - x_bool
            let x_bool_negate = composer.xor_gate(x_bool, one, 10);

            // Define functon part 1: I(x==1)yz
            let function_part_1 = composer.arithmetic_gate(|gate| {
                gate.witness(x_bool, y_times_z, None).mul(F::one())
            });

            // Define functon part 2: I(x!=1)(2y-z)
            let function_part_2 = composer.arithmetic_gate(|gate| {
                gate.witness(x_bool_negate, two_y_minus_z, None)
                    .mul(F::one())
            });

            // Define the full function as
            // r = I(x==1)yz + I(x!=1)(2y-z)
            let full_function = composer.arithmetic_gate(|gate| {
                gate.witness(function_part_1, function_part_2, None)
                    .add(F::one(), F::one())
            });

            composer.assert_equal(full_function, r);

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

    // Prover POV
    let x = 1u64;
    let y = 2u64;
    let z = 5u64;
    let r = 10u64;
    println!("x:{}, y:{}, z:{}, r:{}", x, y, z, r);
    let (proof, pi) = {
        let mut circuit: TestCircuit<BlsScalar, JubJubParameters> =
            TestCircuit {
                x: BlsScalar::from(x),
                y: BlsScalar::from(y),
                z: BlsScalar::from(z),
                r: BlsScalar::from(r),
                dummy: JubJubScalar::from(2u64),
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
