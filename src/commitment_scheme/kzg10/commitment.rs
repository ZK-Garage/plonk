// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Module containing the representation of a Commitment to a Polynomial.
// use std::path::PathBuf;

use ark_ec::PairingEngine;
// use dusk_bytes::{DeserializableSlice, Serializable};

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
/// Holds a commitment to a polynomial in a form of a [`G1Affine`]-bls12_381
/// point.
pub(crate) struct Commitment<E: PairingEngine>(
    /// The commitment is a group element.
    pub(crate) E::G1Affine,
);

impl<E: PairingEngine> From<E::G1Affine> for Commitment<E> {
    fn from(point: E::G1Affine) -> Commitment<E> {
        Commitment(point)
    }
}

impl<E: PairingEngine> From<E::G1Projective> for Commitment<E> {
    fn from(point: E::G1Projective) -> Commitment<E> {
        Commitment(point.into())
    }
}

// impl<E: PairingEngine> Serializable<{ E::G1Affine::SIZE }> for Commitment<E>
// {     type Error = dusk_bytes::Error;

//     fn to_bytes(&self) -> [u8; Self::SIZE] {
//         self.0.to_bytes()
//     }

//     fn from_bytes(buf: &[u8; Self::SIZE]) -> Result<Self, Self::Error> {
//         let g1 = G1Affine::from_slice(buf)?;
//         Ok(Self(g1))
//     }
// }

impl<E: PairingEngine> Commitment<E> {
    /// Builds an identity [`Commitment`] which is equivalent to the
    /// [`G1Affine`] identity point in bls12_381.
    fn identity() -> Commitment<E> {
        Commitment(E::G1Affine::identity())
    }
}

impl<E: PairingEngine> Default for Commitment<E> {
    fn default() -> Commitment<E> {
        Commitment::identity()
    }
}

/*
#[cfg(test)]
mod commitment_tests {
    use super::*;
    use ark_ec::bls12::Bls12;
    use ark_ff::FpParameters;
    use ark_ec::bls12::g1;

    #[test]
    fn commitment_dusk_bytes_serde() {
        let commitment = Commitment(g1::G1Affine);
        let bytes = commitment.to_bytes();
        let obtained_comm = Commitment::from_slice(&bytes)
            .expect("Error on the deserialization");
        assert_eq!(commitment, obtained_comm);
    }
} */
