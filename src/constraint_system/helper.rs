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
use ark_ec::{AffineCurve, PairingEngine, ProjectiveCurve, TEModelParameters};
use ark_ff::{Field, PrimeField, UniformRand};
use ark_poly_commit::kzg10::UniversalParams;
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
    let beta = E::Fr::rand(&mut rng);

    // Compute powers of beta up to and including beta^max_degree
    let powers_of_beta = util::powers_of(&beta, max_degree);

    // Powers of G1 that will be used to commit to a specified polynomial
    let g = E::G1Projective::rand(&mut rng);
    let powers_of_g: Vec<E::G1Affine> =
        slow_multiscalar_mul_single_base(&powers_of_beta, g);
    assert_eq!(powers_of_g.len(), max_degree + 1);

    // Generate the secret scalar gamma
    let gamma = E::Fr::rand(&mut rng);

    // Generate powers of gamma g
    let powers_of_gamma_g_proj: Vec<E::G1Projective> =
        slow_multibase_mul_single_scalar(gamma, &powers_of_g);

    // Compute beta*G2 element and stored cached elements for verifying
    // multiple proofs.
    let h_proj: E::G2Projective = E::G2Projective::rand(&mut rng);
    let h: E::G2Affine = h_proj.into();
    let beta_h: E::G2Affine = (h_proj.mul(beta.into_repr())).into();

    // Compute group elements of the form { \beta^i G2 }
    let neg_powers_of_h: Vec<E::G1Projective> =
        slow_multiscalar_mul_single_base(&powers_of_beta, -h);

    E::G1Projective::batch_normalization(&mut powers_of_gamma_g_proj);
    E::G1Projective::batch_normalization(&mut neg_powers_of_h);

    let powers_of_gamma_g: Vec<E::G1Affine> = powers_of_gamma_g_proj
        .iter()
        .map(|point| point.into_affine())
        .collect();

    let prepared_h: E::G2Prepared = h.into();
    let prepared_beta_h: E::G2Prepared = beta_h.into();

    Ok(UniversalParams {
        powers_of_g: powers_of_g.into(),
        powers_of_gamma_g: powers_of_gamma_g.into(),
        h,
        beta_h,
        neg_powers_of_h: neg_powers_of_h.into(),
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
