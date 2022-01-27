// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) ZK-Garage. All rights reserved.

use crate::lookup::multiset::MultiSet;
use crate::proof_system::widget::{GateConstraint, GateValues};
use ark_ff::{Field, PrimeField};
use core::marker::PhantomData;
use crate::proof_system::linearisation_poly::ProofEvaluations;
use ark_ec::PairingEngine;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::Evaluations;
use ark_poly_commit::sonic_pc::Commitment;
use ark_serialize::*;

/// Lookup Gates Prover Key
// #[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct ProverKey<E> 
where 
E: PairingEngine,
{
    pub(crate) q_lookup: (DensePolynomial<E::Fr>, Evaluations<E::Fr>),
    pub(crate) table_1: (MultiSet<E::Fr>, Commitment<E>, DensePolynomial<E::Fr>),
    pub(crate) table_2: (MultiSet<E::Fr>, Commitment<E>, DensePolynomial<E::Fr>),
    pub(crate) table_3: (MultiSet<E::Fr>, Commitment<E>, DensePolynomial<E::Fr>),
    pub(crate) table_4: (MultiSet<E::Fr>, Commitment<E>, DensePolynomial<E::Fr>),
}

impl<F, E> ProverKey<F, E>
where
E: PairingEngine, 
{

    pub fn compute_quotient_i(&self,
        index: usize,
        w_l_i: E::Fr,
        w_r_i: E::Fr,
        w_o_i: E::Fr,
        w_4_i: E::Fr,
        f_i: F,
        zeta: F,
        lookup_challenge: F,
    ) -> F {


        // q_lookup(X) * (a(X) + zeta * b(X) + (zeta^2 * c(X)) + (zeta^3 * d(X)
        // - f(X))) * α_1
        let a = {
            let q_lookup_i = self.q_lookup.1[index];
            let compressed_tuple =
                Self::compress(w_l_i, w_r_i, w_o_i, w_4_i, zeta);

            q_lookup_i * (compressed_tuple - f_i) * lookup_challenge
        };

        a

    }

    /// Compute linearization for lookup gates
    pub(crate) fn compute_linearization(
        &self,
        a_eval: E::Fr,
        b_eval: E::Fr,
        c_eval: E::Fr,
        d_eval: E::Fr,
        f_eval: E::Fr,
        t_eval: E::Fr,
        t_next_eval: E::Fr,
        h_1_eval: E::Fr,
        h_2_eval: E::Fr,
        p_next_eval: E::Fr,
        l1_eval: E::Fr,
        p_poly: DensePolynomial<E::Fr>,
        h_2_poly: DensePolynomial<E::Fr>,
        (delta, epsilon): (E::Fr, E::Fr),
        zeta: E::Fr,
        lookup_separation_challenge: E::Fr,
    ) -> DensePolynomial<E::Fr> {
        let l_sep_2 = lookup_separation_challenge.square();
        let l_sep_3 = l_sep_2 * lookup_separation_challenge;
        let zeta_sq = zeta * zeta;
        let zeta_cu = zeta * zeta_sq;
        let one_plus_delta = delta + E::Fr::one();
        let epsilon_one_plus_delta = epsilon * one_plus_delta;

        //
        // - q_k(X) * f_eval * lookup_separation_challenge
        let a = {
            let a_0 =
                a_eval + zeta * b_eval + zeta_sq * c_eval + zeta_cu * d_eval;

            &self.q_lookup.0 * ((a_0 - f_eval) * lookup_separation_challenge)
        };

        // p(X) * L0(z) * α_1^2
        let b = { &p_poly * (l1_eval * l_sep_2) };

        // p(X) * (1 + δ) * (ε + f_bar) * (ε(1+δ) + t_bar + δ*tω_bar) * α_1^3
        let c = {
            let c_0 = epsilon + f_eval;
            let c_1 = epsilon_one_plus_delta + t_eval + delta * t_next_eval;

            &p_poly * (one_plus_delta * c_0 * c_1 * l_sep_3)
        };

        // − pω_bar * (ε(1+δ) + h1_bar + δh2_bar) * h2(X) * α_1^3
        let d = {
            let d_0 = epsilon_one_plus_delta + h_1_eval + delta * h_2_eval;

            &h_2_poly * (-p_next_eval * d_0 * l_sep_3)
        };

        a + b + c + d
    }

    fn compress(
        w_l: E::Fr,
        w_r: E::Fr,
        w_o: E::Fr,
        w_4: E::Fr,
        zeta: E::Fr,
    ) -> E::Fr {
        let zeta_sq = zeta.square();
        let zeta_cu = zeta_sq * zeta;
    
        let a = w_l;
    
        let b = w_r * zeta;
    
        let c = w_o * zeta_sq;
    
        let d = w_4 * zeta_cu;
    
        a + b + c + d
    }
    
}

