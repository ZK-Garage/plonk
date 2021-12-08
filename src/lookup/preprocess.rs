// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Module containing the lookup functionality

use crate::error::Error;
use crate::lookup::{LookupTable, MultiSet};
use ark_ec::PairingEngine;
use ark_poly::domain::EvaluationDomain;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly_commit::kzg10::{Commitment, Powers, KZG10};

/// This table will be the preprocessed version of the
/// precomputed table, T, with arity 4. This structure
/// is passed to the proof alongside the table of witness
/// values.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct PreprocessedLookupTable<E>
where
    E: PairingEngine,
{
    /// This is the circuit size
    pub n: u32,

    /// This is the first column in the preprocessed
    /// table containing a MultiSet, Commitments to the
    /// MultiSet and the coefficients as a Polynomial
    pub(crate) t_1: (MultiSet<E::Fr>, Commitment<E>, DensePolynomial<E::Fr>),

    /// This is the second column in the preprocessed
    /// table containing a MultiSet, Commitments to the
    /// MultiSet and the coefficients as a Polynomial
    pub(crate) t_2: (MultiSet<E::Fr>, Commitment<E>, DensePolynomial<E::Fr>),

    /// This is the third column in the preprocessed
    /// table containing a MultiSet, Commitments to the
    /// MultiSet and the coefficients as a Polynomial
    pub(crate) t_3: (MultiSet<E::Fr>, Commitment<E>, DensePolynomial<E::Fr>),

    /// This is the fourth column in the preprocessed
    /// table containing a MultiSet, Commitments to the
    /// MultiSet and the coefficients as a Polynomial
    pub(crate) t_4: (MultiSet<E::Fr>, Commitment<E>, DensePolynomial<E::Fr>),
}

impl<E> PreprocessedLookupTable<E>
where
    E: PairingEngine,
{
    /// This function takes in a precomputed look up table and
    /// pads it to the length of the circuit entries, as a power
    /// of 2. The function then interpolates a polynomial from the
    /// padded table and makes a commitment to the poly. The
    /// outputted struct will be used in the proof alongside our
    /// circuit witness table.
    pub fn preprocess(
        table: &LookupTable<E::Fr>,
        commit_key: &Powers<E>,
        n: u32,
    ) -> Result<Self, Error> {
        let domain = EvaluationDomain::new(n as usize).unwrap();

        let columned_table = table.vec_to_multiset();
        let mut t_1 = columned_table.0;
        let mut t_2 = columned_table.1;
        let mut t_3 = columned_table.2;
        let mut t_4 = columned_table.3;

        t_1.pad(n);
        t_2.pad(n);
        t_3.pad(n);
        t_4.pad(n);

        let t_1_poly = t_1.to_polynomial(&domain);
        let t_2_poly = t_2.to_polynomial(&domain);
        let t_3_poly = t_3.to_polynomial(&domain);
        let t_4_poly = t_4.to_polynomial(&domain);

        let t_1_commit = KZG10::commit(commit_key, &t_1_poly, None, None)?.0;
        let t_2_commit = KZG10::commit(commit_key, &t_2_poly, None, None)?.0;
        let t_3_commit = KZG10::commit(commit_key, &t_3_poly, None, None)?.0;
        let t_4_commit = KZG10::commit(commit_key, &t_4_poly, None, None)?.0;

        Ok(Self {
            n,
            t_1: (t_1, t_1_commit, t_1_poly),
            t_2: (t_2, t_2_commit, t_2_poly),
            t_3: (t_3, t_3_commit, t_3_poly),
            t_4: (t_4, t_4_commit, t_4_poly),
        })
    }
}
