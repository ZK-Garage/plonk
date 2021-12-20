// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use ark_ff::Field;
use crate::error::Error;
//use crate::fft::{EvaluationDomain, Polynomial};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
    UVPolynomial, Polynomial,
};
use core::ops::{Add, Mul};
use ark_ec::PairingEngine;


/// MultiSet is struct containing vectors of scalars, which
/// individually represents either a wire value or an index
/// of a PlookUp table
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct MultiSet<E: PairingEngine>(pub Vec<E::Fr>);

impl<E: PairingEngine> Default for MultiSet<E> {
    fn default() -> MultiSet<E> {
        MultiSet::new()
    }
}

impl<E: PairingEngine> From<&[E::Fr]> for MultiSet<E> {
    fn from(slice: &[E::Fr]) -> MultiSet<E> {
        MultiSet(slice.to_vec())
    }
}

impl<E: PairingEngine> 
    MultiSet<E> {
    /// Creates an empty vector with a multiset wrapper around it
    pub fn new() -> MultiSet<E> {
        MultiSet(vec![])
        
    }

    /// Generate a `MultiSet` struct from a slice of bytes.
    pub fn from_slice(bytes: &[u8]) -> Result<MultiSet<E>, Error> {
        let mut buffer = bytes;
        let elements = buffer
            .chunks(E::Fr::SIZE)
            .map(|chunk| E::Fr::from_slice(chunk))
            .collect::<Result<Vec<E::Fr>, Error>>()?;
        Ok(MultiSet(elements))
    }

    /// Given a [`MultiSet`], return it in it's bytes representation
    /// element by element.
    pub fn to_var_bytes(&self) -> Vec<u8> {
        self.0
            .iter()
            .map(|item| item.to_bytes().to_vec())
            .flatten()
            .collect()
    }

    /// Extends the length of the multiset to n elements
    /// The n will be the size of the arithmetic circuit.
    /// This will extend the vectors to the size
    pub fn pad(&mut self, n: u32) {
        assert!(n.is_power_of_two());
        let diff = n - self.len() as u32;
        self.0.extend(vec![self.0[0]; diff as usize]);
    }

    /// Pushes chosen value onto the end of the Multiset
    pub fn push(&mut self, value: E::Fr) {
        self.0.push(value)
    }

    /// Fetches last element in MultiSet.
    /// Returns None if there are no elements in the MultiSet.
    pub fn last(&self) -> Option<&E::Fr> {
        self.0.last()
    }

    /// Returns the cardinality of the multiset
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Returns whether or not the multiset is empty.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// Returns the position of the element in the Multiset.
    /// Returns None if the element is not found.
    pub fn position(&self, element: &E::Fr) -> Option<usize> {
        self.0.iter().position(|&x| x == *element)
    }

    /// Concatenates and sorts two Multisets together.
    /// From the Plookup paper, if we have t: {1,2,4,3}
    /// and f: {2,3,4,1}.
    /// We first check if all elements of f exist in t
    /// Then we combine the multisets together and sort
    /// their elements together. The final MultiSet will
    /// look as follows, s: {1,1,2,2,3,3,4,4}
    pub fn sorted_concat(&self, f: &MultiSet<E>) -> Result<MultiSet<E>, Error> {
        let mut s = self.clone();
        s.0.reserve(f.0.len());
        for element in f.0.iter() {
            let index = s.position(element).ok_or(Error::ElementNotIndexed)?;
            s.0.insert(index, *element);
        }

        Ok(s)
    }

    /// Checks whether one mutltiset is a subset of another.
    /// This function will be used to check if the all elements
    /// in set f, from the paper, are contained inside t.
    pub fn contains_all(&self, other: &MultiSet<E>) -> bool {
        other.0.iter().all(|item| self.contains(item))
    }

    /// Checks if an element is in the MultiSet
    pub fn contains(&self, entry: &E::Fr) -> bool {
        self.0.contains(entry)
    }

    /// Splits a multiset into halves as specified by the paper
    /// The last element of the first half should be the same
    /// as the first element of the second half.
    /// Since a multiset can never have an even cardinality, we
    /// always split it in the way described above.
    pub fn halve(&self) -> (MultiSet<E>, MultiSet<E>) {
        let length = self.0.len();

        let first_half = MultiSet::from(&self.0[0..=length / 2]);
        let second_half = MultiSet::from(&self.0[length / 2..]);

        (first_half, second_half)
    }

    /// Splits a multiset into alternating halves of the same length
    /// as specified in the Plonkup paper. A multiset must have even
    /// cardinality to be split in this manner.
    pub fn halve_alternating(&self) -> (MultiSet<E>, MultiSet<E>) {
        let mut evens = vec![];
        let mut odds = vec![];
        for i in 0..self.len() {
            if i % 2 == 0 {
                evens.push(self.0[i]);
            } else {
                odds.push(self.0[i]);
            }
        }

        (MultiSet(evens), MultiSet(odds))
    }

    /// Treats each element in the multiset as evaluation points
    /// Computes IFFT of the set of evaluation points
    /// and returns the coefficients as a Polynomial data structure
    pub(crate) fn to_polynomial(
        &self,
        domain: &GeneralEvaluationDomain<E::Fr>,
    ) -> DensePolynomial<E::Fr> {
        DensePolynomial::from_coefficients_vec(domain.ifft(&self.0))
    }

    /// Turn three multisets into a single multiset using
    /// a random challenge, Alpha. Alpha is dervived by hashing
    /// the transcript.
    /// The function iterates over the given sets and mutiplies by alpha:
    /// a + (b * alpha) + (c * alpha^2)  
    pub fn compress_three_arity(
        multisets: [&MultiSet<E>; 3],
        alpha: E::Fr,
    ) -> MultiSet<E> {
        MultiSet(
            multisets[0]
                .0
                .iter()
                .zip(multisets[1].0.iter())
                .zip(multisets[2].0.iter())
                .map(|((a, b), c)| a + b * alpha + c * alpha.square())
                .collect::<Vec<E::Fr>>(),
        )
    }

    /// Turn four multisets into a single multiset using
    /// a random challenge, Alpha. Alpha is dervived by hashing
    /// the transcript.
    /// The function iterates over the given sets and mutiplies by alpha:
    /// a + (b * alpha) + (c * alpha^2) + (d * alpha^3)  
    pub fn compress_four_arity(
        multisets: [&MultiSet<E>; 4],
        alpha: E::Fr,
    ) -> MultiSet<E> {
        MultiSet(
            multisets[0]
                .0
                .iter()
                .zip(multisets[1].0.iter())
                .zip(multisets[2].0.iter())
                .zip(multisets[3].0.iter())
                .map(|(((a, b), c), d)| {
                    a + b * alpha
                        + c * alpha.square()
                        + d * alpha.pow(&[3u64, 0u64, 0u64, 0u64])
                })
                .collect::<Vec<E::Fr>>(),
            )
    }
}

impl<E: PairingEngine> Add for MultiSet<E> {
    type Output = MultiSet<E>;

    fn add(self, other: MultiSet<E>) -> Self::Output {
        let result = self
            .0
            .into_iter()
            .zip(other.0.iter())
            .map(|(x, y)| x + y)
            .collect();

        MultiSet(result)
    }
}

impl<E: PairingEngine> Mul for MultiSet<E> {
    type Output = MultiSet<E>;

    fn mul(self, other: MultiSet<E>) -> Self::Output {
        let result = self
            .0
            .into_iter()
            .zip(other.0.iter())
            .map(|(x, y)| x * y)
            .collect();

        MultiSet(result)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_poly::EvaluationDomain;
    use crate::lookup::WitnessTable;

    #[test]
    fn test_halve<E: PairingEngine>() {
        let mut s = MultiSet::new();
        s.push(E::Fr::from(0));
        s.push(E::Fr::from(1));
        s.push(E::Fr::from(2));
        s.push(E::Fr::from(3));
        s.push(E::Fr::from(4));
        s.push(E::Fr::from(5));
        s.push(E::Fr::from(6));

        let (h_1, h_2) = s.halve();
        assert_eq!(h_1.len(), 4);
        assert_eq!(h_2.len(), 4);

        let left_half = MultiSet(vec![
            E::Fr::from(0),
            E::Fr::from(1),
            E::Fr::from(2),
            E::Fr::from(3),
        ]);

        assert_eq!(left_half, h_1);

        let right_half = MultiSet(vec![
            E::Fr::from(3),
            E::Fr::from(4),
            E::Fr::from(5),
            E::Fr::from(6),
        ]);

        assert_eq!(right_half, h_2);

        // The last element of the first half should equal the first
        // element of the second half.
        assert_eq!(h_1.0.last().unwrap(), &h_2.0[0])
    }

    #[test]
    fn test_to_polynomial<E: PairingEngine>() {
        let mut s = MultiSet::new();
        s.push(E::Fr::from(1));
        s.push(E::Fr::from(2));
        s.push(E::Fr::from(3));
        s.push(E::Fr::from(4));
        s.push(E::Fr::from(5));
        s.push(E::Fr::from(6));
        s.push(E::Fr::from(7));

        let domain = EvaluationDomain::new(s.len() + 1).unwrap();
        let s_poly = s.to_polynomial(&domain);

        assert_eq!(s_poly.degree(), 7)
    }
    #[test]
    fn test_is_subset<E: PairingEngine>() {
        let mut t = MultiSet::new();
        t.push(E::Fr::from(1));
        t.push(E::Fr::from(2));
        t.push(E::Fr::from(3));
        t.push(E::Fr::from(4));
        t.push(E::Fr::from(5));
        t.push(E::Fr::from(6));
        t.push(E::Fr::from(7));
        let mut f = MultiSet::new();
        f.push(E::Fr::from(1));
        f.push(E::Fr::from(2));
        let mut n = MultiSet::new();
        n.push(E::Fr::from(8));

        assert!(t.contains_all(&f));
        assert!(!t.contains_all(&n));
    }

    #[test]
    fn test_full_compression_into_s<E: PairingEngine>() {
        let mut t = MultiSet::new();

        t.push(E::Fr::zero());
        t.push(E::Fr::one());
        t.push(E::Fr::from(2));
        t.push(E::Fr::from(3));
        t.push(E::Fr::from(4));
        t.push(E::Fr::from(5));
        t.push(E::Fr::from(6));
        t.push(E::Fr::from(7));

        let mut f = MultiSet::new();
        f.push(E::Fr::from(3));
        f.push(E::Fr::from(6));
        f.push(E::Fr::from(0));
        f.push(E::Fr::from(5));
        f.push(E::Fr::from(4));
        f.push(E::Fr::from(3));
        f.push(E::Fr::from(2));
        f.push(E::Fr::from(0));
        f.push(E::Fr::from(0));
        f.push(E::Fr::from(1));
        f.push(E::Fr::from(2));

        assert!(t.contains_all(&f));

        assert!(t.contains(&E::Fr::from(2)));

        let s = t.sorted_concat(&f);

        // The sets should be merged but also
        // in the ascending order
        let concatenated_set = MultiSet(vec![
            E::Fr::zero(),
            E::Fr::zero(),
            E::Fr::zero(),
            E::Fr::zero(),
            E::Fr::one(),
            E::Fr::one(),
            E::Fr::from(2),
            E::Fr::from(2),
            E::Fr::from(2),
            E::Fr::from(3),
            E::Fr::from(3),
            E::Fr::from(3),
            E::Fr::from(4),
            E::Fr::from(4),
            E::Fr::from(5),
            E::Fr::from(5),
            E::Fr::from(6),
            E::Fr::from(6),
            E::Fr::from(7),
        ]);

        assert_eq!(s.unwrap(), concatenated_set);
    }

    #[test]
    fn multiset_compression_input<E: PairingEngine>() {
        // Alpha is a random challenge from
        // the transcript
        let alpha = E::Fr::from(2);
        let alpha_squared = alpha * alpha;

        let mut table = WitnessTable::default();

        // Fill in wires directly, no need to use a
        // plookup table as this will not be going
        // into a proof
        table.from_wire_values(
            E::Fr::from(1),
            E::Fr::from(2),
            E::Fr::from(3),
            E::Fr::from(3),
        );

        // Computed expected result
        let compressed_element = MultiSet::compress_three_arity(
            [&table.f_1, &table.f_2, &table.f_3],
            alpha,
        );

        let actual_element = E::Fr::from(1)
            + (E::Fr::from(2) * alpha)
            + (E::Fr::from(3) * alpha_squared);

        let mut actual_set = MultiSet::new();

        actual_set.push(actual_element);

        assert_eq!(actual_set, compressed_element);
    }
}
