// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::constraint_system::StandardComposer;
use crate::error::Error;
use crate::proof_system::widget::VerifierKey as PlonkVerifierKey;
use crate::proof_system::Proof;
use crate::transcript::TranscriptWrapper;
use ark_ec::{PairingEngine, TEModelParameters};
use ark_poly_commit::kzg10::{Powers, VerifierKey};

/// Abstraction structure designed verify [`Proof`]s.
#[allow(missing_debug_implementations)]
pub struct Verifier<
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
> {
    /// VerificationKey which is used to verify a specific PLONK circuit
    pub verifier_key: Option<PlonkVerifierKey<E, P>>,

    pub(crate) cs: StandardComposer<E, P>,
    /// Store the messages exchanged during the preprocessing stage
    /// This is copied each time, we make a proof, so that we can use the same
    /// verifier to Verify multiple proofs from the same circuit. If this
    /// is not copied, then the verification procedure will modify
    /// the transcript, making it unusable for future proofs.
    pub preprocessed_transcript: TranscriptWrapper<E>,
}

impl<
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    > Default for Verifier<E, P>
{
    fn default() -> Verifier<E, P> {
        Verifier::new(b"plonk")
    }
}

impl<
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    > Verifier<E, P>
{
    /// Creates a new `Verifier` instance.
    pub fn new(label: &'static [u8]) -> Verifier<E, P> {
        Verifier {
            verifier_key: None,
            cs: StandardComposer::new(),
            preprocessed_transcript: TranscriptWrapper::new(label),
        }
    }

    /// Creates a new `Verifier` instance with some expected size.
    pub fn with_expected_size(
        label: &'static [u8],
        size: usize,
    ) -> Verifier<E, P> {
        Verifier {
            verifier_key: None,
            cs: StandardComposer::with_expected_size(size),
            preprocessed_transcript: TranscriptWrapper::new(label),
        }
    }

    /// Returns the number of gates in the circuit.
    pub fn circuit_size(&self) -> usize {
        self.cs.circuit_size()
    }

    /// Returns a mutable copy of the underlying composer.
    pub fn mut_cs(&mut self) -> &mut StandardComposer<E, P> {
        &mut self.cs
    }

    /// Preprocess a circuit to obtain a [`PlonkVerifierKey`] and a circuit
    /// descriptor so that the `Verifier` instance can verify [`Proof`]s
    /// for this circuit descriptor instance.
    pub fn preprocess(&mut self, commit_key: &Powers<E>) -> Result<(), Error> {
        let vk = self.cs.preprocess_verifier(
            commit_key,
            &mut self.preprocessed_transcript,
        )?;

        self.verifier_key = Some(vk);
        Ok(())
    }

    /// Keys the [`Transcript`] with additional seed information
    /// Wrapper around [`Transcript::append_message`].
    pub fn key_transcript(&mut self, label: &'static [u8], message: &[u8]) {
        self.preprocessed_transcript
            .transcript
            .append_message(label, message);
    }

    /// Verifies a [`Proof`].
    pub fn verify(
        &self,
        proof: &Proof<E, P>,
        pc_verifier_key: &VerifierKey<E>,
        public_inputs: &[E::Fr],
    ) -> Result<(), Error> {
        let mut cloned_transcript = self.preprocessed_transcript.clone();
        let plonk_verifier_key = self.verifier_key.as_ref().unwrap();

        proof.verify(
            plonk_verifier_key,
            &mut cloned_transcript,
            pc_verifier_key,
            public_inputs,
        )
    }
}
