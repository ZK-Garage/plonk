// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! PLONK Example for stateless interface

use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
use ark_ec::{TEModelParameters, models::twisted_edwards_extended::GroupAffine, AffineCurve, ProjectiveCurve};
use ark_ff::PrimeField;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::{sonic_pc::SonicKZG10, PolynomialCommitment};
use rand::rngs::OsRng;
use ark_ed_on_bls12_381::{
    EdwardsParameters as JubJubParameters, Fr as JubJubScalar,
};

use plonk::{
    constraint_system::StandardComposer, 
    prelude::{Error, from_embedded_curve_scalar}, 
    proof_system::interface::{compile, prove, verify}, circuit::PublicInputBuilder,
};

// Implements a circuit that checks:
// 1) a + b = c where C is a PI
// 2) a <= 2^6
// 3) b <= 2^5
// 4) a * b = d where D is a PI
// 5) JubJub::GENERATOR * e(JubJubScalar) = f where F is a PI
fn generate_circuit<F, P>(
    composer: &mut StandardComposer<F, P>,
    a: F,
    b: F,
    c: F,
    d: F,
    e: P::ScalarField,
    f: GroupAffine<P>,
) -> Result<(), Error> 
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    let a = composer.add_input(a);
    let b = composer.add_input(b);
    let zero = composer.zero_var();

    // Make first constraint a + b = c (as public input)
    composer.arithmetic_gate(|gate| {
        gate.witness(a, b, Some(zero))
            .add(F::one(), F::one())
            .pi(-c)
    });

    // Check that a and b are in range
    composer.range_gate(a, 1 << 6);
    composer.range_gate(b, 1 << 5);
    // Make second constraint a * b = d
    composer.arithmetic_gate(|gate| {
        gate.witness(a, b, Some(zero)).mul(F::one()).pi(-d)
    });
    let e =
        composer.add_input(from_embedded_curve_scalar::<F, P>(e));
    let (x, y) = P::AFFINE_GENERATOR_COEFFS;
    let generator = GroupAffine::new(x, y);
    let scalar_mul_result =
        composer.fixed_base_scalar_mul(e, generator);
    // Apply the constrain
    composer.assert_equal_public_point(scalar_mul_result, f);
    Ok(())
}

fn main() -> Result<(), Error> {
    type PC = SonicKZG10<Bls12_381, DensePolynomial<BlsScalar>>;

    // Generate circuit
    let circuit_size = 1 << 11;    
    let mut composer = StandardComposer::<BlsScalar, JubJubParameters>::default();
    let (x, y) = JubJubParameters::AFFINE_GENERATOR_COEFFS;
    let generator: GroupAffine<JubJubParameters> = GroupAffine::new(x, y);
    let point_f_pi: GroupAffine<JubJubParameters> =
        AffineCurve::mul(&generator, JubJubScalar::from(2u64).into_repr())
            .into_affine();
    generate_circuit::<BlsScalar, JubJubParameters>(
        &mut composer,
        BlsScalar::from(20u64),
        BlsScalar::from(5u64),
        BlsScalar::from(25u64),
        BlsScalar::from(100u64),
        JubJubScalar::from(2u64),
        point_f_pi,
    ).expect("Can't generate circuit");

    // Compile
    let pp = PC::setup(1 << 12, None, &mut OsRng)
        .expect("Unable to sample public parameters.");

    let (proving_key, verifying_key) = compile::<BlsScalar, JubJubParameters, PC>(
        &pp,
        &mut composer,
        circuit_size,
    ).unwrap();

    // Prove
    let proof = prove::<BlsScalar, JubJubParameters, PC>(
        &proving_key,
        composer,
        circuit_size,
    ).unwrap();

    // Verify
    let public_inputs = PublicInputBuilder::new()
        .add_input(&BlsScalar::from(25u64))
        .unwrap()
        .add_input(&BlsScalar::from(100u64))
        .unwrap()
        .add_input(&point_f_pi)
        .unwrap()
        .finish();

    let verify_status = verify::<BlsScalar, JubJubParameters, PC>(
        &verifying_key,
        &public_inputs,
        &proof,
    ).unwrap();
    println!("Successfully verify: {}", verify_status);

    Ok(())
}
