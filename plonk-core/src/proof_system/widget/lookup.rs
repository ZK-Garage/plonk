// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) ZK-Garage. All rights reserved.

use crate::lookup::multiset::MultiSet;
use crate::proof_system::widget::{GateConstraint, GateValues};
use ark_ff::{Field, PrimeField, FftField};
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
pub struct ProverKey<F> 
where 
F: PrimeField,
{
    pub q_lookup: (DensePolynomial<F>, Evaluations<F>),
    pub table_1: (MultiSet<F>, Commitment<F>, DensePolynomial<F>),
    pub table_2: (MultiSet<F>, Commitment<F>, DensePolynomial<F>),
    pub table_3: (MultiSet<F>, Commitment<F>, DensePolynomial<F>),
    pub table_4: (MultiSet<F>, Commitment<F>, DensePolynomial<F>),
}

impl<F> ProverKey<F>
where
F: PrimeField, 
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
        zeta: F,
        lookup_separation_challenge: F,
    ) -> DensePolynomial<F> {
        let l_sep_2 = lookup_separation_challenge.square();
        let l_sep_3 = l_sep_2 * lookup_separation_challenge;
        let zeta_sq = zeta * zeta;
        let zeta_cu = zeta * zeta_sq;

        // q_lookup(X) * f_eval * lookup_separation_challenge
        let a = {
            let a_0 =
                a_eval + zeta * b_eval + zeta_sq * c_eval + zeta_cu * d_eval;

            &self.q_lookup.0 * ((a_0 - f_eval) * l_sep_3)
        };

        a

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

