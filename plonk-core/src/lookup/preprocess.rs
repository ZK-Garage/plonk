// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::commitment::HomomorphicCommitment;
use crate::error::{Error, to_pc_error};
use crate::label_polynomial;
use crate::lookup::{LookupTable, MultiSet};
use ark_ec::PairingEngine;
use ark_ff::PrimeField;
use ark_poly::domain::EvaluationDomain;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly_commit::kzg10::{Commitment, Powers, KZG10};

/// This table will be the preprocessed version of the
/// precomputed table, T, with arity 4. This structure
/// is passed to the proof alongside the table of witness
/// values.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct PreprocessedLookupTable<F, PC>
where
F: PrimeField,
PC: HomomorphicCommitment<F>,
{
    /// This is the circuit size
    pub n: u32,

    /// This is the first column in the preprocessed
    /// table containing a MultiSet, Commitments to the
    /// MultiSet and the coefficients as a Polynomial
    pub(crate) t_1: (MultiSet<F>, PC::Commitment, DensePolynomial<F>),

    /// This is the second column in the preprocessed
    /// table containing a MultiSet, Commitments to the
    /// MultiSet and the coefficients as a Polynomial
    pub(crate) t_2: (MultiSet<F>, PC::Commitment, DensePolynomial<F>),

    /// This is the third column in the preprocessed
    /// table containing a MultiSet, Commitments to the
    /// MultiSet and the coefficients as a Polynomial
    pub(crate) t_3: (MultiSet<F>, PC::Commitment, DensePolynomial<F>),

    /// This is the fourth column in the preprocessed
    /// table containing a MultiSet, Commitments to the
    /// MultiSet and the coefficients as a Polynomial
    pub(crate) t_4: (MultiSet<F>, PC::Commitment, DensePolynomial<F>),
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
        // n should be a power of 2
        assert!(n & (n - 1) == 0);

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

        let (table_column_commitments, _) = PC::commit(commit_key, &[label_polynomial!(t_1_poly), label_polynomial!(t_2_poly), label_polynomial!(t_3_poly), label_polynomial!(t_4_poly), ] , None).map_err(to_pc_error::<F, PC>)?;

        Ok(Self {
            n,
            t_1: (t_1, table_column_commitments[0].commitment().clone(), t_1_poly),
            t_2: (t_2, table_column_commitments[1].commitment().clone(), t_2_poly),
            t_3: (t_3, table_column_commitments[2].commitment().clone(), t_3_poly),
            t_4: (t_4, table_column_commitments[3].commitment().clone(), t_4_poly),
        })
    }
}

#[cfg(test)]
mod test {
    use crate::batch_test;
    use crate::commitment::HomomorphicCommitment;
    use crate::constraint_system::StandardComposer;
    use crate::error::to_pc_error;
    use crate::lookup::{LookupTable, PreprocessedLookupTable};
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_ec::{PairingEngine, TEModelParameters};
    use ark_ff::PrimeField;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly_commit::{
        kzg10::{Powers, KZG10},
        sonic_pc::SonicKZG10,
        PolynomialCommitment,
    };
    use rand::rngs::OsRng;

    /// This function creates a table and preprocesses it. Then it checks that
    /// all table columns are the same length.
    fn test_table_preprocessing<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let universal_params =
        PC::setup(32, None, &mut OsRng).map_err(to_pc_error::<F, PC>).unwrap();

        // Commit Key
        let (ck, _) = PC::trim(
            &universal_params,
            32,
            0,
            None,
        )
        .map_err(to_pc_error::<F, PC>).unwrap();

        let mut table: LookupTable<F> = LookupTable::new();

        (0..11).for_each(|a| {
            table.insert_xor_row(19u64, 6u64, 64u64);
            table.insert_xor_row(4u64, 2u64, 64u64);
        });

        let preprocessed_table =
            PreprocessedLookupTable::<F, PC>::preprocess(&table, &ck, 32).unwrap();

        assert!(
            preprocessed_table.n as usize == preprocessed_table.t_1.0.len()
        );
        assert!(
            preprocessed_table.n as usize == preprocessed_table.t_2.0.len()
        );
        assert!(
            preprocessed_table.n as usize == preprocessed_table.t_3.0.len()
        );
        assert!(
            preprocessed_table.n as usize == preprocessed_table.t_4.0.len()
        );
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
