// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::proof_system::linearisation_poly::ProofEvaluations;
use crate::proof_system::widget::range::proverkey::delta;
use ark_ec::PairingEngine;
use ark_ff::Field;
use ark_poly_commit::sonic_pc::Commitment;
use ark_serialize::*;
#[derive(
    Debug, PartialEq, Eq, Copy, Clone, CanonicalDeserialize, CanonicalSerialize,
)]
pub(crate) struct VerifierKey<E: PairingEngine> {
    pub(crate) q_range: Commitment<E>,
}

impl<E: PairingEngine> VerifierKey<E> {
    pub(crate) fn compute_linearisation_commitment(
        &self,
        range_separation_challenge: E::Fr,
        scalars: &mut Vec<E::Fr>,
        points: &mut Vec<E::G1Affine>,
        evaluations: &ProofEvaluations<E::Fr>,
    ) {
        let four = E::Fr::from(4u64);

        let kappa = range_separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;

        let b_1 = delta(evaluations.c_eval - (four * evaluations.d_eval));
        let b_2 = delta(evaluations.b_eval - four * evaluations.c_eval) * kappa;
        let b_3 =
            delta(evaluations.a_eval - four * evaluations.b_eval) * kappa_sq;
        let b_4 = delta(evaluations.d_next_eval - (four * evaluations.a_eval))
            * kappa_cu;

        scalars.push((b_1 + b_2 + b_3 + b_4) * range_separation_challenge);
        points.push(self.q_range.0);
    }
}
