// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) ZK-Garage. All rights reserved.

//! Public Inputs of the circuit. This values are available for the
//! [`super::Prover`] and [`super::Verifier`].
//!
//! This module contains the implementation of the [`PI`] struct and all the
//! basic manipulations such as inserting new values and getting the public
//! inputs in evaluation or coefficient form.

use core::ops::Deref;

use alloc::collections::BTreeMap;
use ark_ff::{PrimeField, ToConstraintField};
use ark_poly::{
    polynomial::UVPolynomial, univariate::DensePolynomial, EvaluationDomain,
    GeneralEvaluationDomain,
};
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write,
};

use crate::prelude::Error;

///  Public Inputs
#[derive(
    CanonicalDeserialize, CanonicalSerialize, Clone, PartialEq, Debug, Hash,
)]
pub struct PI<F>
where
    F: PrimeField,
{
    // padded size of the circuit
    n: usize,
    // non-zero values of the public input
    values: BTreeMap<usize, F>,
}

impl<F> Deref for PI<F>
where
    F: PrimeField,
{
    type Target = BTreeMap<usize, F>;
    fn deref(&self) -> &Self::Target {
        &self.values
    }
}

impl<F> PI<F>
where
    F: PrimeField,
{
    /// Creates a new struct for [`PI`].
    pub fn new(n: usize) -> Self {
        assert!(n.is_power_of_two());
        Self {
            n,
            values: BTreeMap::new(),
        }
    }

    /// Updates the size of the circuit.
    pub fn update_size(&mut self, n: usize) {
        assert!(n.is_power_of_two());
        self.n = n;
    }

    /// Inserts a new public input value at a given position.
    ///
    /// If the value provided is zero no insertion will be made as zeros are the
    /// implicit value of empty positions.
    /// The function will panic if an insertion is atempted in an already
    /// occupied position.
    pub fn insert(&mut self, pos: usize, val: F) {
        if self.values.contains_key(&pos) {
            panic!("Insertion in public inputs conflicts with previous value at position {}", pos);
        }
        if val != F::zero() {
            self.values.insert(pos, val);
        }
    }

    /// Inserts public input data that can be converted to one or more field
    /// elements starting at a given position.
    /// Returns the number of field elements occupied by the input or
    /// [`Error::InvalidPublicInputValue`] if the input could not be converted.
    pub fn add_input<T>(&mut self, pos: usize, item: &T) -> Result<usize, Error>
    where
        T: ToConstraintField<F>,
    {
        self.extend(pos, [item])
    }

    /// Inserts all the elements of `iter` converted into constraint field
    /// elements in consecutive positions.
    /// Returns the number of field elements occupied by `iter`
    /// [`Error::InvalidPublicInputValue`] if the input could not be
    /// converted.
    fn extend<'t, T, I>(
        &mut self,
        init_pos: usize,
        iter: I,
    ) -> Result<usize, Error>
    where
        T: 't + ToConstraintField<F>,
        I: IntoIterator<Item = &'t T>,
    {
        let mut result = 0;
        for (pos1, item) in iter.into_iter().enumerate() {
            let item_repr = &item
                .to_field_elements()
                .ok_or(Error::InvalidPublicInputValue)?;

            for (pos2, elem) in item_repr.iter().enumerate() {
                self.insert(init_pos + pos1 + pos2, *elem);
                result += 1;
            }
        }
        Ok(result)
    }

    /// Returns the public inputs as a vector of `n` evaluations.
    pub fn as_evals(&self) -> Vec<F> {
        let mut pi = vec![F::zero(); self.n];
        self.values.iter().for_each(|(pos, eval)| pi[*pos] = *eval);
        pi
    }
}

impl<F> From<&PI<F>> for DensePolynomial<F>
where
    F: PrimeField,
{
    fn from(pi: &PI<F>) -> Self {
        let domain = GeneralEvaluationDomain::<F>::new(pi.n).unwrap();
        let evals = pi.as_evals();
        DensePolynomial::from_coefficients_vec(domain.ifft(&evals))
    }
}
