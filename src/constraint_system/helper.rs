// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use super::StandardComposer;
use crate::error::Error;
use crate::proof_system::{Prover, Verifier};
use ark_ec::{PairingEngine, ProjectiveCurve, TEModelParameters};
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::kzg10::{self, Powers, KZG10};
use ark_poly_commit::sonic_pc::SonicKZG10;
use ark_poly_commit::PolynomialCommitment;
use num_traits::{One, Zero};
use rand_core::OsRng;

/// Adds dummy constraints using arithmetic gates
pub(crate) fn dummy_gadget<
    E: PairingEngine,
    T: ProjectiveCurve<BaseField = E::Fr>,
    P: TEModelParameters<BaseField = E::Fr>,
>(
    n: usize,
    composer: &mut StandardComposer<E, T, P>,
) {
    let one = E::Fr::one();

    let var_one = composer.add_input(one);

    for _ in 0..n {
        composer.big_add(
            (E::Fr::one(), var_one),
            (E::Fr::one(), var_one),
            None,
            E::Fr::zero(),
            None,
        );
    }
}

/// Takes a generic gadget function with no auxillary input and
/// tests whether it passes an end-to-end test
pub(crate) fn gadget_tester<
    E: PairingEngine,
    T: ProjectiveCurve<BaseField = E::Fr>,
    P: TEModelParameters<BaseField = E::Fr>,
>(
    gadget: fn(composer: &mut StandardComposer<E, T, P>),
    n: usize,
) -> Result<(), Error> {
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
        gadget(&mut prover.mut_cs());

        // Commit Key
        let (ck, _) = SonicKZG10::<E, DensePolynomial<E::Fr>>::trim(
            &universal_params,
            prover.cs.circuit_size().next_power_of_two(),
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
    gadget(&mut verifier.mut_cs());

    // Compute Commit and Verifier Key
    let (sonic_ck, sonic_vk) = SonicKZG10::<E, DensePolynomial<E::Fr>>::trim(
        &universal_params,
        verifier.cs.circuit_size().next_power_of_two(),
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
