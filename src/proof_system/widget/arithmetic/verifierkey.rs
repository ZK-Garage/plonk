// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use ark_ec::PairingEngine;
use ark_poly_commit::sonic_pc::Commitment;
use ark_serialize::*;
#[derive(
    Debug, PartialEq, Eq, Copy, Clone, CanonicalDeserialize, CanonicalSerialize,
)]
pub(crate) struct VerifierKey<E: PairingEngine> {
    pub q_m: Commitment<E>,
    pub q_l: Commitment<E>,
    pub q_r: Commitment<E>,
    pub q_o: Commitment<E>,
    pub q_4: Commitment<E>,
    pub q_c: Commitment<E>,
    pub q_arith: Commitment<E>,
}

use crate::proof_system::linearisation_poly::ProofEvaluations;

impl<E: PairingEngine> VerifierKey<E> {
    pub(crate) fn compute_linearisation_commitment(
        &self,
        scalars: &mut Vec<E::Fr>,
        points: &mut Vec<E::G1Affine>,
        evaluations: &ProofEvaluations<E::Fr>,
    ) {
        let q_arith_eval = evaluations.q_arith_eval;

        scalars.push(evaluations.a_eval * &evaluations.b_eval * &q_arith_eval);
        points.push(self.q_m.0);

        scalars.push(evaluations.a_eval * &q_arith_eval);
        points.push(self.q_l.0);

        scalars.push(evaluations.b_eval * &q_arith_eval);
        points.push(self.q_r.0);

        scalars.push(evaluations.c_eval * &q_arith_eval);
        points.push(self.q_o.0);

        scalars.push(evaluations.d_eval * &q_arith_eval);
        points.push(self.q_4.0);

        scalars.push(q_arith_eval);
        points.push(self.q_c.0);
    }
}
