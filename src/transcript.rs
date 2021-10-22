// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! This is an extension over the [Merlin Transcript](Transcript)
//! which adds a few extra functionalities.
use ark_ec::PairingEngine;
use ark_ff::{Field, PrimeField};
use ark_poly_commit::kzg10::Commitment;
use ark_serialize::CanonicalSerialize;
use core::marker::PhantomData;
use merlin::Transcript;

#[derive(Clone)]
pub struct TranscriptWrapper<E: PairingEngine> {
    pub transcript: Transcript,
    _marker: PhantomData<E>,
}

impl<E: PairingEngine> TranscriptWrapper<E> {
    pub fn new(label: &'static [u8]) -> TranscriptWrapper<E> {
        TranscriptWrapper {
            transcript: Transcript::new(label),
            _marker: PhantomData,
        }
    }
}

/// Transcript adds an abstraction over the Merlin transcript
/// For convenience
pub(crate) trait TranscriptProtocol<E: PairingEngine> {
    /// Append a `commitment` with the given `label`.
    fn append_commitment(&mut self, label: &'static [u8], comm: &Commitment<E>);

    /// Append a scalar with the given `label`.
    fn append_scalar(&mut self, label: &'static [u8], s: &E::Fr);

    /// Compute a `label`ed challenge variable.
    fn challenge_scalar(&mut self, label: &'static [u8]) -> E::Fr;

    /// Append domain separator for the circuit size.
    fn circuit_domain_sep(&mut self, n: u64);
}

impl<E: PairingEngine> TranscriptProtocol<E> for TranscriptWrapper<E> {
    fn append_commitment(
        &mut self,
        label: &'static [u8],
        comm: &Commitment<E>,
    ) {
        let mut bytes = Vec::new();
        comm.0.serialize(&mut bytes).unwrap();
        self.transcript.append_message(label, &bytes);
    }

    fn append_scalar(&mut self, label: &'static [u8], s: &E::Fr) {
        let mut bytes = Vec::new();
        s.serialize(&mut bytes).unwrap();
        self.transcript.append_message(label, &bytes)
    }

    fn challenge_scalar(&mut self, label: &'static [u8]) -> E::Fr {
        // XXX: review this
        let size = E::Fr::size_in_bits() / 8;
        let mut buf = vec![0u8; size];
        self.transcript.challenge_bytes(label, &mut buf);
        // XXX: This may fail
        E::Fr::from_random_bytes(&buf).unwrap()
    }

    fn circuit_domain_sep(&mut self, n: u64) {
        self.transcript.append_message(b"dom-sep", b"circuit_size");
        self.transcript.append_u64(b"n", n);
    }
}
