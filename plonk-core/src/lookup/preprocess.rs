// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::commitment::HomomorphicCommitment;
use crate::error::{to_pc_error, Error};
use crate::lookup::{LookupTable, MultiSet};
use ark_ff::PrimeField;
use ark_poly::domain::EvaluationDomain;
use ark_poly::polynomial::univariate::DensePolynomial;

/// This table will be the preprocessed version of the precomputed table,
/// T, with arity 4. This structure is passed to the proof alongside the
/// table of witness values.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PreprocessedLookupTable<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    /// This is the circuit size
    pub n: u32,

    /// Vector of columns in the preprocessed table containing a
    /// `MultiSet`, `Commitment` and `DensePolynomial`.
    pub(crate) t: Vec<(MultiSet<F>, PC::Commitment, DensePolynomial<F>)>,
}

impl<F, PC> PreprocessedLookupTable<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    /// This function takes in a precomputed look up table and
    /// pads it to the length of the circuit entries, as a power
    /// of 2. The function then interpolates a polynomial from the
    /// padded table and makes a commitment to the poly. The
    /// outputted struct will be used in the proof alongside our
    /// circuit witness table.
    pub fn preprocess(
        table: &LookupTable<F>,
        commit_key: &PC::CommitterKey,
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
                let commitment = PC::commit(commit_key, &[labeled_poly], None)
                    .map_err(to_pc_error::<F, PC>)?;
                Ok((column, commitment.0[0].commitment().clone(), poly))
            })
            .collect::<Result<_, Error>>()?;
        Ok(Self { n, t: result })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::batch_test;
    use crate::commitment::HomomorphicCommitment;
    use crate::lookup::{LookupTable, PreprocessedLookupTable};
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_ec::TEModelParameters;
    use rand_core::OsRng;

    /// This function creates a table and preprocesses it. Then it checks that
    /// all table columns are the same length.
    #[allow(clippy::extra_unused_type_parameters)]
    fn test_table_preprocessing<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let pp = PC::setup(32, None, &mut OsRng)
            .map_err(to_pc_error::<F, PC>)
            .unwrap();
        let (_committer_key, _) = PC::trim(&pp, 32, 0, None)
            .map_err(to_pc_error::<F, PC>)
            .unwrap();

        // Commit Key
        let (ck, _) = PC::trim(&pp, 32, 0, None)
            .map_err(to_pc_error::<F, PC>)
            .unwrap();

        let mut table: LookupTable<F> = LookupTable::new();

        (0..11).for_each(|_a| {
            table.insert_xor_row(19u64, 6u64, 64u64);
            table.insert_xor_row(4u64, 2u64, 64u64);
        });

        let preprocessed_table =
            PreprocessedLookupTable::<F, PC>::preprocess(&table, &ck, 32)
                .unwrap();

        preprocessed_table.t.iter().for_each(|column| {
            assert!(preprocessed_table.n as usize == column.0.len());
        });
    }

    // Bls12-381 tests
    batch_test!(
        [
            test_table_preprocessing
        ],
        [] => (
            Bls12_381,
            ark_ed_on_bls12_381::EdwardsParameters
        )
    );

    // Bls12-377 tests
    batch_test!(
        [
            test_table_preprocessing
        ],
        [] => (
            Bls12_377,
            ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
