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
pub struct ProverKey<F, E> 
where 
F: PrimeField, 
E: PairingEngine,
{
    pub(crate) q_lookup: (DensePolynomial<F>, Evaluations<F>),
    pub(crate) table_1: (MultiSet<E::Fr>, Commitment<E>, DensePolynomial<F>),
    pub(crate) table_2: (MultiSet<E::Fr>, Commitment<E>, DensePolynomial<F>),
    pub(crate) table_3: (MultiSet<E::Fr>, Commitment<E>, DensePolynomial<F>),
    pub(crate) table_4: (MultiSet<E::Fr>, Commitment<E>, DensePolynomial<F>),
}

impl<F, E> ProverKey<F, E>
where
F: PrimeField, 
E: PairingEngine, 
{

    pub fn compute_quotient_i(&self,
        index: usize,
        w_l_i: F,
        w_r_i: F,
        w_o_i: F,
        w_4_i: F,
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
        a_eval: F,
        b_eval: F,
        c_eval: F,
        d_eval: F,
        f_eval: F,
        t_eval: F,
        t_next_eval: F,
        h_1_eval: F,
        h_2_eval: F,
        p_next_eval: F,
        l1_eval: F,
        p_poly: DensePolynomial<F>,
        h_2_poly: DensePolynomial<F>,
        (delta, epsilon): (F, F),
        zeta: F,
        lookup_separation_challenge: F,
    ) -> DensePolynomial<F> {
        let l_sep_2 = lookup_separation_challenge.square();
        let l_sep_3 = l_sep_2 * lookup_separation_challenge;
        let zeta_sq = zeta * zeta;
        let zeta_cu = zeta * zeta_sq;
        let one_plus_delta = delta + F::one();
        let epsilon_one_plus_delta = epsilon * one_plus_delta;

        //
        // - q_k(X) * f_eval * lookup_separation_challenge
        let a = {
            let a_0 =
                a_eval + zeta * b_eval + zeta_sq * c_eval + zeta_cu * d_eval;

            &self.q_lookup.0 * &((a_0 - f_eval) * lookup_separation_challenge)
        };

        // p(X) * L0(z) * α_1^2
        let b = { p_poly * &(l1_eval * l_sep_2) };

        // p(X) * (1 + δ) * (ε + f_bar) * (ε(1+δ) + t_bar + δ*tω_bar) * α_1^3
        let c = {
            let c_0 = epsilon + f_eval;
            let c_1 = epsilon_one_plus_delta + t_eval + delta * t_next_eval;

            p_poly * &(one_plus_delta * c_0 * c_1 * l_sep_3)
        };

        // − pω_bar * (ε(1+δ) + h1_bar + δh2_bar) * h2(X) * α_1^3
        let d = {
            let d_0 = epsilon_one_plus_delta + h_1_eval + delta * h_2_eval;

            &(-p_next_eval * d_0 * l_sep_3) * h_2_poly;
        };

        let mut r = a;
        r += &b;
        r += &c;
        r += &d;

        r
    }

    fn compress(
        w_l: F,
        w_r: F,
        w_o: F,
        w_4: F,
        zeta: F,
    ) -> F {
        let zeta_sq = zeta.square();
        let zeta_cu = zeta_sq * zeta;
    
        let a = w_l;
    
        let b = w_r * zeta;
    
        let c = w_o * zeta_sq;
    
        let d = w_4 * zeta_cu;
    
        a + b + c + d
    }
    
}

