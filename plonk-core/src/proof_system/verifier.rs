// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Verifier-side of the PLONK Proving System

//use crate::circuit::EmbeddedCurve;
use super::pi::PublicInputs;
use crate::{
    constraint_system::StandardComposer,
    error::Error,
    parameters::CircuitParameters,
    proof_system::{widget::VerifierKey as PlonkVerifierKey, Proof},
};
use merlin::Transcript;

/// Abstraction structure designed verify [`Proof`]s.
pub struct Verifier<P>
where
    P: CircuitParameters,
{
    /// VerificationKey which is used to verify a specific PLONK circuit
    pub verifier_key: Option<PlonkVerifierKey<P>>,

    /// Circuit Description
    pub(crate) cs: StandardComposer<P>,

    /// Store the messages exchanged during the preprocessing stage.
    ///
    /// This is copied each time, we make a proof, so that we can use the same
    /// verifier to verify multiple proofs from the same circuit. If this is
    /// not copied, then the verification procedure will modify the transcript,
    /// making it unusable for future proofs.
    pub preprocessed_transcript: Transcript,
}

impl<P> Verifier<P>
where
    P: CircuitParameters,
{
    /// Creates a new `Verifier` instance.
    pub fn new(label: &'static [u8]) -> Self {
        Self {
            verifier_key: None,
            cs: StandardComposer::new(),
            preprocessed_transcript: Transcript::new(label),
        }
    }

    /// Creates a new `Verifier` instance with some expected size.
    pub fn with_expected_size(label: &'static [u8], size: usize) -> Self {
        Self {
            verifier_key: None,
            cs: StandardComposer::with_expected_size(size),
            preprocessed_transcript: Transcript::new(label),
        }
    }

    /// Returns the smallest power of two needed for the curcuit
    pub fn circuit_bound(&self) -> usize {
        self.cs.circuit_bound()
    }

    /// Returns a mutable copy of the underlying composer.
    pub fn mut_cs(&mut self) -> &mut StandardComposer<P> {
        &mut self.cs
    }

    /// Preprocess a circuit to obtain a [`PlonkVerifierKey<F, PC>`] and a
    /// circuit descriptor so that the `Verifier` instance can verify
    /// [`Proof`]s for this circuit descriptor instance.
    pub fn preprocess(
        &mut self,
        commit_key: &P::CommitterKey,
    ) -> Result<(), Error> {
        let vk = self.cs.preprocess_verifier(
            commit_key,
            &mut self.preprocessed_transcript,
        )?;

        self.verifier_key = Some(vk);
        Ok(())
    }

    /// Keys the [`Transcript`] with additional seed information
    /// Wrapper around [`Transcript::append_message`].
    ///
    /// [`Transcript`]: merlin::Transcript
    /// [`Transcript::append_message`]: merlin::Transcript::append_message
    pub fn key_transcript(&mut self, label: &'static [u8], message: &[u8]) {
        self.preprocessed_transcript.append_message(label, message);
    }

    /// Verifies a [`Proof`] using `pc_verifier_key` and `public_inputs`.
    pub fn verify(
        &self,
        proof: &Proof<P>,
        pc_verifier_key: &P::VerifierKey,
        public_inputs: &PublicInputs<P::ScalarField>,
    ) -> Result<(), Error> {
        proof.verify(
            self.verifier_key.as_ref().unwrap(),
            &mut self.preprocessed_transcript.clone(),
            pc_verifier_key,
            public_inputs,
        )
    }
}

impl<P> Default for Verifier<P>
where
    P: CircuitParameters,
{
    #[inline]
    fn default() -> Verifier<P> {
        Verifier::new(b"plonk")
    }
}
