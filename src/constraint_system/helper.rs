// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use super::StandardComposer;
use crate::error::Error;
use crate::proof_system::{Prover, Verifier};
use ark_ec::{ModelParameters, PairingEngine, TEModelParameters};
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::kzg10::{self, Powers, KZG10};
use ark_poly_commit::sonic_pc::SonicKZG10;
use ark_poly_commit::PolynomialCommitment;
use rand_core::OsRng;

use ark_ff::FftField;

/// Adds dummy constraints using arithmetic gates.
#[allow(dead_code)]
pub(crate) fn dummy_gadget<F, P>(
    n: usize,
    composer: &mut StandardComposer<F, P>,
) where
    F: FftField,
    P: ModelParameters<BaseField = F>,
{
    let one = F::one();
    let var_one = composer.add_input(one);
    for _ in 0..n {
        composer.big_add(
            (F::one(), var_one),
            (F::one(), var_one),
            None,
            F::zero(),
            None,
        );
    }
}

/// Takes a generic gadget function with no auxillary input and tests whether it
/// passes an end-to-end test.
#[allow(dead_code)]
pub(crate) fn gadget_tester<E, P>(
    gadget: fn(&mut StandardComposer<E::Fr, P>),
    n: usize,
) -> Result<(), Error>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    // Common View
    let universal_params =
        SonicKZG10::<E, DensePolynomial<E::Fr>>::setup(2 * n, None, &mut OsRng)
            .unwrap();

    type PC<E> = SonicKZG10<E, DensePolynomial<<E as PairingEngine>::Fr>>;
    // Provers View
    let (proof, public_inputs) = {
        // Create a prover struct
        let mut prover = Prover::<E::Fr, P, PC<E>>::new(b"demo");

        // Additionally key the transcript
        prover.key_transcript(b"key", b"additional seed information");

        // Add gadgets
        gadget(prover.mut_cs());

        // Commit Key
        let (ck, _) = PC::<E>::trim(
            &universal_params,
            prover.circuit_size().next_power_of_two(),
            0,
            None,
        )
        .unwrap();

        // Preprocess circuit
        prover.preprocess(&ck)?;

        // Once the prove method is called, the public inputs are cleared
        // So pre-fetch these before calling Prove
        let public_inputs = prover.cs.construct_dense_pi_vec();

        // Compute Proof
        (prover.prove(&ck)?, public_inputs)
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
    let (ck, vk) = PC::<E>::trim(
        &universal_params,
        verifier.circuit_size().next_power_of_two(),
        0,
        None,
    )
    .unwrap();

    // Preprocess circuit
    verifier.preprocess(&ck)?;

    // Verify proof
    verifier.verify(&proof, &vk, &public_inputs)
}
