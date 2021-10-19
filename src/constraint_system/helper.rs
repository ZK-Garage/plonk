// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use super::StandardComposer;
// use ark_poly_commit::PublicParameters;
use crate::error::Error;
use crate::proof_system::{Prover, Verifier};
use crate::util;
use ark_ec::{PairingEngine, ProjectiveCurve, TEModelParameters};
use ark_poly_commit::kzg10::UniversalParams;
use rand_core::{CryptoRng, OsRng, RngCore};

/// Setup generates the public parameters using a random number generator.
/// This method will in most cases be used for testing and exploration.
/// In reality, a `Trusted party` or a `Multiparty Computation` will be used
/// to generate the SRS. Returns an error if the configured degree is less
/// than one.
pub fn setup<R: RngCore + CryptoRng, E: PairingEngine>(
    max_degree: usize,
    mut rng: &mut R,
) -> Result<UniversalParams<E>, Error> {
    // Cannot commit to constants
    if max_degree < 1 {
        return Err(Error::DegreeIsZero);
    }

    // Generate the secret scalar beta
    let beta = util::random_scalar(&mut rng);

    // Compute powers of beta up to and including beta^max_degree
    let powers_of_beta = util::powers_of(&beta, max_degree);

    // Powers of G1 that will be used to commit to a specified polynomial
    let g = util::random_g1_point(&mut rng);
    let powers_of_g: Vec<E::G1Projective> =
        util::slow_multiscalar_mul_single_base(&powers_of_beta, g);
    assert_eq!(powers_of_g.len(), max_degree + 1);

    // Generate the secret scalar gamma
    let gamma = util::random_scalar(&mut rng);

    // Generate powers of gamma g
    let powers_of_gamma_g: Vec<E::G1Projective> =
        util::slow_multibase_mul_single_scalar(gamma, &powers_of_g);

    // Normalise all projective points
    let mut normalised_g = vec![E::G1Affine::identity(); max_degree + 1];
    E::G1Projective::batch_normalize(&powers_of_g, &mut normalised_g);

    let mut normalised_gamma_g = vec![E::G1Affine::identity(); max_degree + 1];
    E::G1Projective::batch_normalize(
        &powers_of_gamma_g,
        &mut normalised_gamma_g,
    );

    // Compute beta*G2 element and stored cached elements for verifying
    // multiple proofs.
    let h: E::G2Affine = util::random_g2_point(&mut rng).into();
    let beta_h: E::G2Affine = (h * beta).into();

    // Compute group elements of the form { \beta^i G2 }
    let powers_of_h: Vec<E::G1Projective> =
        util::slow_multiscalar_mul_single_base(&powers_of_beta, h);

    let mut normalised_h = vec![E::G1Affine::identity(); max_degree + 1];
    E::G1Projective::batch_normalize(&powers_of_h, &mut normalised_h);

    let prepared_h: E::G2Prepared = h.into();
    let prepared_beta_h: E::G2Prepared = beta_h.into();

    Ok(UniversalParams {
        powers_of_g: normalised_g,
        powers_of_gamma_g: normalised_gamma_g,
        h,
        beta_h,
        neg_powers_of_h: powers_of_h,
        prepared_h,
        prepared_beta_h,
    })
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
    let unviersal_params = UniversalParams::setup(2 * n, &mut OsRng)?;
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
