// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use super::StandardComposer;
// use ark_poly_commit::PublicParameters;
use crate::commitment_scheme::kzg10::PublicParameters;
use crate::error::Error;
use crate::proof_system::{Prover, Verifier};
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve, TEModelParameters};
use ark_ff::{Field, PrimeField, UniformRand};
use num_traits::{One, Zero};
use rand_core::{CryptoRng, OsRng, RngCore};

/// This function is only used to generate the SRS.
/// The intention is just to compute the resulting points
/// of the operation `a*P, b*P, c*P ... (n-1)*P` into a `Vec`.
pub(crate) fn slow_multiscalar_mul_single_base<E: PairingEngine>(
    scalars: &[E::Fr],
    base: E::G1Projective,
) -> Vec<E::G1Projective> {
    scalars.iter().map(|&s| base.mul(s.into_repr())).collect()
}

/// This function is only used to generate the SRS.
pub(crate) fn slow_multibase_mul_single_scalar<E: PairingEngine>(
    scalar: E::Fr,
    bases: &[E::G1Projective],
) -> Vec<E::G1Projective> {
    let scalar_repr = scalar.into_repr();
    bases.iter().map(|&b| b.mul(scalar_repr)).collect()
}

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
    let unviersal_params = PublicParameters::setup(2 * n, &mut OsRng)?;
    // Provers View
    let (proof, public_inputs) = {
        // Create a prover struct
        let mut prover = Prover::new(b"demo");

        // Additionally key the transcript
        prover.key_transcript(b"key", b"additional seed information");

        // Add gadgets
        gadget(&mut prover.mut_cs());

        // Commit Key
        let (ck, _) = unviersal_params
            .trim(2 * prover.cs.circuit_size().next_power_of_two())?;

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
    gadget(&mut verifier.mut_cs());

    // Compute Commit and Verifier Key
    let (ck, vk) = unviersal_params
        .trim(verifier.cs.circuit_size().next_power_of_two())?;

    // Preprocess circuit
    verifier.preprocess(&ck)?;

    // Verify proof
    verifier.verify(&proof, &vk, &public_inputs)
}
