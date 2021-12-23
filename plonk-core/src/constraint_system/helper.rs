// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! Testing Helper Functions

use super::StandardComposer;
use crate::error::Error;
use crate::proof_system::{Prover, Verifier};
use ark_ec::{PairingEngine, TEModelParameters};
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::kzg10::{self, Powers, KZG10};
use ark_poly_commit::sonic_pc::SonicKZG10;
use ark_poly_commit::PolynomialCommitment;
use num_traits::One;
use rand_core::OsRng;

/// Adds dummy constraints using arithmetic gates.
#[allow(dead_code)]
pub(crate) fn dummy_gadget<E, P>(
    n: usize,
    composer: &mut StandardComposer<E, P>,
) where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    let one = E::Fr::one();
    let var_one = composer.add_input(one);
    for _ in 0..n {
        composer.arithmetic_gate(|gate| {
            gate.witness(var_one, var_one, None)
                .add(E::Fr::one(), E::Fr::one())
        });
    }
}

/// Takes a generic gadget function with no auxillary input and tests whether it
/// passes an end-to-end test.
#[allow(dead_code)]
pub(crate) fn gadget_tester<E, P>(
    gadget: fn(&mut StandardComposer<E, P>),
    n: usize,
) -> Result<(), Error>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    // Common View
    let universal_params =
        KZG10::<E, DensePolynomial<E::Fr>>::setup(2 * n, false, &mut OsRng)?;
    // Provers View
    let (proof, public_inputs) = {
        // Create a prover struct
        let mut prover = Prover::new(b"demo");

        // Additionally key the transcript
        prover.key_transcript(b"key", b"additional seed information");

        // Add gadgets
        gadget(prover.mut_cs());

        // Commit Key
        let (ck, _) = SonicKZG10::<E, DensePolynomial<E::Fr>>::trim(
            &universal_params,
            prover.circuit_size().next_power_of_two(),
            0,
            None,
        )
        .unwrap();
        let powers = Powers {
            powers_of_g: ck.powers_of_g.into(),
            powers_of_gamma_g: ck.powers_of_gamma_g.into(),
        };
        // Preprocess circuit
        prover.preprocess(&powers)?;

        // Once the prove method is called, the public inputs are cleared
        // So pre-fetch these before calling Prove
        let public_inputs = prover.cs.construct_dense_pi_vec();

        // Compute Proof
        (prover.prove(&powers)?, public_inputs)
    };
    // Verifiers view
    //
    // Create a Verifier object
    let mut verifier = Verifier::new(b"demo");

    // Additionally key the transcript
    verifier.key_transcript(b"key", b"additional seed information");

    // Add gadgets
    gadget(verifier.mut_cs());

    // Compute Commit and Verifier Key
    let (sonic_ck, sonic_vk) = SonicKZG10::<E, DensePolynomial<E::Fr>>::trim(
        &universal_params,
        verifier.circuit_size().next_power_of_two(),
        0,
        None,
    )
    .unwrap();
    let powers = Powers {
        powers_of_g: sonic_ck.powers_of_g.into(),
        powers_of_gamma_g: sonic_ck.powers_of_gamma_g.into(),
    };

    let vk = kzg10::VerifierKey {
        g: sonic_vk.g,
        gamma_g: sonic_vk.gamma_g,
        h: sonic_vk.h,
        beta_h: sonic_vk.beta_h,
        prepared_h: sonic_vk.prepared_h,
        prepared_beta_h: sonic_vk.prepared_beta_h,
    };
    // Preprocess circuit
    verifier.preprocess(&powers)?;

    // Verify proof
    verifier.verify(&proof, &vk, &public_inputs)
}
