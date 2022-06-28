// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::{
    error::{to_pc_error, Error},
    lookup::{LookupTable, MultiSet},
    parameters::CircuitParameters,
};
use ark_poly::{
    domain::EvaluationDomain, polynomial::univariate::DensePolynomial,
};
use ark_poly_commit::PolynomialCommitment;

/// This table will be the preprocessed version of the precomputed table,
/// T, with arity 4. This structure is passed to the proof alongside the
/// table of witness values.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PreprocessedLookupTable<P>
where
    P: CircuitParameters,
{
    /// This is the circuit size
    pub n: u32,

    /// Vector of columns in the preprocessed table containing a
    /// `MultiSet`, `Commitment` and `DensePolynomial`.
    pub(crate) t: Vec<(
        MultiSet<P::ScalarField>,
        P::Commitment,
        DensePolynomial<P::ScalarField>,
    )>,
}

impl<P> PreprocessedLookupTable<P>
where
    P: CircuitParameters,
{
    /// This function takes in a precomputed look up table and
    /// pads it to the length of the circuit entries, as a power
    /// of 2. The function then interpolates a polynomial from the
    /// padded table and makes a commitment to the poly. The
    /// outputted struct will be used in the proof alongside our
    /// circuit witness table.
    pub fn preprocess(
        table: &LookupTable<P::ScalarField>,
        commit_key: &P::CommitterKey,
        n: u32,
    ) -> Result<Self, Error> {
        assert!(n.is_power_of_two());
        let domain = EvaluationDomain::new(n as usize).unwrap();
        let result = table
            .vec_to_multiset()
            .into_iter()
            .enumerate()
            .map(|(index, mut column)| {
                column.pad(n);
                let poly = column.to_polynomial(&domain);
                let poly_name = format!("t_{}_poly", index + 1);
                let labeled_poly = ark_poly_commit::LabeledPolynomial::new(
                    poly_name,
                    poly.to_owned(),
                    None,
                    None,
                );
                let commitment = P::PolynomialCommitment::commit(
                    commit_key,
                    &[labeled_poly],
                    None,
                )
                .map_err(to_pc_error::<P>)?;
                Ok((column, commitment.0[0].commitment().clone(), poly))
            })
            .collect::<Result<_, Error>>()?;
        Ok(Self { n, t: result })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        batch_test,
        lookup::{LookupTable, PreprocessedLookupTable},
        parameters::test::*,
    };
    use rand_core::OsRng;

    /// This function creates a table and preprocesses it. Then it checks that
    /// all table columns are the same length.
    fn test_table_preprocessing<P>()
    where
        P: CircuitParameters,
    {
        let pp = P::PolynomialCommitment::setup(32, None, &mut OsRng)
            .map_err(to_pc_error::<P>)
            .unwrap();
        let (_committer_key, _) =
            P::PolynomialCommitment::trim(&pp, 32, 0, None)
                .map_err(to_pc_error::<P>)
                .unwrap();

        // Commit Key
        let (ck, _) = P::PolynomialCommitment::trim(&pp, 32, 0, None)
            .map_err(to_pc_error::<P>)
            .unwrap();

        let mut table: LookupTable<P::ScalarField> = LookupTable::new();

        (0..11).for_each(|_a| {
            table.insert_xor_row(19u64, 6u64, 64u64);
            table.insert_xor_row(4u64, 2u64, 64u64);
        });

        let preprocessed_table =
            PreprocessedLookupTable::<P>::preprocess(&table, &ck, 32).unwrap();

        preprocessed_table.t.iter().for_each(|column| {
            assert!(preprocessed_table.n as usize == column.0.len());
        });
    }

    // Bls12-381 tests
    batch_test!(
        [
            test_table_preprocessing
        ],
        [] => [
            Bls12_381_KZG, Bls12_381_IPA, Bls12_377_KZG, Bls12_377_IPA
        ]
    );
}
