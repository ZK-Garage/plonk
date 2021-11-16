// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::permutation::constants::{K1, K2, K3};
use crate::proof_system::linearisation_poly::ProofEvaluations;
use ark_ec::PairingEngine;
use ark_ff::Field;
use ark_poly_commit::sonic_pc::Commitment;
use ark_serialize::*;
#[derive(
    Debug, PartialEq, Eq, Copy, Clone, CanonicalDeserialize, CanonicalSerialize,
)]
pub(crate) struct VerifierKey<E: PairingEngine> {
    pub(crate) left_sigma: Commitment<E>,
    pub(crate) right_sigma: Commitment<E>,
    pub(crate) out_sigma: Commitment<E>,
    pub(crate) fourth_sigma: Commitment<E>,
}

impl<E: PairingEngine> VerifierKey<E> {
    pub(crate) fn new(
        left_sigma: Commitment<E>,
        right_sigma: Commitment<E>,
        out_sigma: Commitment<E>,
        fourth_sigma: Commitment<E>,
    ) -> VerifierKey<E> {
        VerifierKey {
            left_sigma,
            right_sigma,
            out_sigma,
            fourth_sigma,
        }
    }
    pub(crate) fn compute_linearisation_commitment(
        &self,
        scalars: &mut Vec<E::Fr>,
        points: &mut Vec<E::G1Affine>,
        evaluations: &ProofEvaluations<E::Fr>,
        z_challenge: E::Fr,
        (alpha, beta, gamma): (E::Fr, E::Fr, E::Fr),
        l1_eval: E::Fr,
        z_comm: E::G1Affine,
    ) {
        let alpha_sq = alpha.square();

        // (a_eval + beta * z + gamma)(b_eval + beta * z * k1 +
        // gamma)(c_eval + beta * k2 * z + gamma)(d_eval + beta
        // * k3 * z + gamma) * alpha
        let x = {
            let beta_z = beta * z_challenge;
            let q_0 = evaluations.a_eval + beta_z + gamma;

            let beta_k1_z = beta * K1::<E::Fr>() * z_challenge;
            let q_1 = evaluations.b_eval + beta_k1_z + gamma;

            let beta_k2_z = beta * K2::<E::Fr>() * z_challenge;
            let q_2 = evaluations.c_eval + beta_k2_z + gamma;

            let beta_k3_z = beta * K3::<E::Fr>() * z_challenge;
            let q_3 = (evaluations.d_eval + beta_k3_z + gamma) * alpha;

            q_0 * q_1 * q_2 * q_3
        };

        // l1(z) * alpha^2
        let r = l1_eval * alpha_sq;

        scalars.push(x + r);
        points.push(z_comm);

        // -(a_eval + beta * sigma_1_eval + gamma)(b_eval + beta *
        // sigma_2_eval + gamma)(c_eval + beta * sigma_3_eval +
        // gamma) * alpha^2
        let y = {
            let beta_sigma_1 = beta * evaluations.left_sigma_eval;
            let q_0 = evaluations.a_eval + beta_sigma_1 + gamma;

            let beta_sigma_2 = beta * evaluations.right_sigma_eval;
            let q_1 = evaluations.b_eval + beta_sigma_2 + gamma;

            let beta_sigma_3 = beta * evaluations.out_sigma_eval;
            let q_2 = evaluations.c_eval + beta_sigma_3 + gamma;

            let q_3 = beta * evaluations.perm_eval * alpha;

            -(q_0 * q_1 * q_2 * q_3)
        };
        scalars.push(y);
        points.push(self.fourth_sigma.0);
    }
}
