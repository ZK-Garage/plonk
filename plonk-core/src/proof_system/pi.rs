// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) ZK-Garage. All rights reserved.

//! Public Inputs of the circuit. This values are available for the `Prover` and
//! `Verifier`.
//!
//! This module contains the implementation of the `PI` struct and all the basic
//! manipulations such as inserting new values and getting the public inputs in
//! evaluation or coefficient form.

use ark_ff::PrimeField;
use ark_poly::{
    polynomial::UVPolynomial, univariate::DensePolynomial, EvaluationDomain,
    GeneralEvaluationDomain,
};

///  Public Inputs
#[derive(Clone, PartialEq)]
pub struct PI<F>
where
    F: PrimeField,
{
    // padded size of the circuit
    n: usize,
    // non-zero values of the public input
    values: Vec<(usize, F)>,
}

impl<F> PI<F>
where
    F: PrimeField,
{
    //TODO Remove if unused
    // #[inline]
    // pub fn len(&self) -> usize {
    //     self.values.len()
    // }

    // #[inline]
    // pub fn is_empty(&self) -> bool {
    //     self.values.is_empty()
    // }

    // #[inline]
    // pub fn iter<'a>(&'a self) -> Iter<'a, F> {
    //     Iter::new(&self.values[..])
    // }

    /// Creates a new struct for `PI`.
    pub fn new(n: usize) -> Self {
        assert!(n.is_power_of_two());
        Self { n, values: vec![] }
    }

    /// Inserts a new public input value.
    ///
    /// Only the last value for any given position is taken into account.
    /// Inserting a value effectively overrides any other value inserted at
    /// the same position.
    pub fn insert(&mut self, pos: usize, val: F) {
        assert!(pos < self.n);
        self.values.push((pos, val));
    }

    /// Returns the public inputs as a vector of `n` evaluations.
    pub fn as_evals(&self) -> Vec<F> {
        let mut pi = vec![F::zero(); self.n];
        self.values.iter().for_each(|(pos, eval)| pi[*pos] = *eval);
        pi
    }
}

impl<F> Into<DensePolynomial<F>> for PI<F>
where
    F: PrimeField,
{
    fn into(self) -> DensePolynomial<F> {
        let domain = GeneralEvaluationDomain::<F>::new(self.n).unwrap();
        let evals = self.as_evals();
        DensePolynomial::from_coefficients_vec(domain.fft(&evals))
    }
}
