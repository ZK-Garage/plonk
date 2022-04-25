// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::error::Error;
use crate::lookup::{LookupTable, MultiSet};
use ark_ff::Field;

/// This witness table contains queries to a lookup table for lookup gates.
/// This table can have any arity. (But for the time other parts of the codebase
/// force 4-arity)
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct WitnessTable<F>
where
    F: Field,
{
    /// Vector containing one `MultiSet` for each wire
    pub f: Vec<MultiSet<F>>,
}

impl<F> WitnessTable<F>
where
    F: Field,
{
    /// Initializes empty witness table of arity 4
    pub fn new() -> Self {
        Self {
            f: vec![MultiSet::new(); 4],
        }
    }

    /// This allows the witness table to be filled directly without
    /// taking any values, or the results, from the lookup table.
    /// If the values do no exist in the lookup table, then the proof
    /// will fail when witness and preprocessed tables are concatenated.
    pub fn from_wire_values(&mut self, wires: Vec<F>) {
        assert_eq!(wires.len(), self.f.len());
        wires
            .into_iter()
            .zip(self.f.iter_mut())
            .for_each(|(val, column)| column.push(val));
    }

    /// Attempts to look up a value from a lookup table. If successful, all four
    /// elements are pushed to their respective multisets.
    pub fn value_from_table(
        &mut self,
        lookup_table: &LookupTable<F>,
        left_wire_val: F,
        right_wire_val: F,
        fourth_wire_val: F,
    ) -> Result<(), Error> {
        // This function is specific for 4-arity
        assert_eq!(self.f.len(), 4);
        let output_wire_val = lookup_table.lookup(
            left_wire_val,
            right_wire_val,
            fourth_wire_val,
        )?;
        self.f[0].push(left_wire_val);
        self.f[1].push(right_wire_val);
        self.f[2].push(output_wire_val);
        self.f[3].push(fourth_wire_val);
        Ok(())
    }
}
#[cfg(test)]
mod test {
    use super::*;
    use crate::batch_field_test;
    use crate::lookup::LookupTable;
    use ark_bls12_377::Fr as bls12_377_scalar_field;
    use ark_bls12_381::Fr as bls12_381_scalar_field;

    fn test_lookup_fuctionality_1<F>()
    where
        F: Field,
    {
        // Build lookup table
        let lookup_table = LookupTable::<F>::xor_table(0, 3);

        // Instantiate empty multisets of wire values in witness table
        let mut f = WitnessTable::<F>::new();
        // Read values from lookup table and insert into witness table
        assert!(f
            .value_from_table(
                &lookup_table,
                F::from(2u64),
                F::from(5u64),
                -F::one()
            )
            .is_ok());
        // Check that non existent elements cause a failure
        assert!(f
            .value_from_table(
                &lookup_table,
                F::from(25u64),
                F::from(5u64),
                -F::one()
            )
            .is_err());
    }

    fn test_lookup_fuctionality_2<F>()
    where
        F: Field,
    {
        // Build empty lookup tables
        let mut lookup_table = LookupTable::<F>::new();

        // Add a consecutive set of tables, with
        // XOR operationd and addition operations
        lookup_table.insert_multi_xor(0, 4);
        lookup_table.insert_multi_add(2, 3);

        // Build empty witness table
        let mut f = WitnessTable::<F>::new();

        // Check for output of wires within lookup table and
        // if they exist input them to the witness table
        assert!(f
            .value_from_table(
                &lookup_table,
                F::from(2u32),
                F::from(3u32),
                -F::one()
            )
            .is_ok());
        assert!(f
            .value_from_table(
                &lookup_table,
                F::from(4u32),
                F::from(6u32),
                F::zero()
            )
            .is_ok());

        // Check that values not contained in the lookup table
        // do not get added to the witness table
        assert!(f
            .value_from_table(
                &lookup_table,
                F::from(22u32),
                F::one(),
                -F::one()
            )
            .is_err());
        assert!(f
            .value_from_table(&lookup_table, F::zero(), F::one(), F::zero())
            .is_err());
    }

    // Bls12-381 tests
    batch_field_test!(
        [
            test_lookup_fuctionality_1,
            test_lookup_fuctionality_2
        ],
        [] => bls12_381_scalar_field
    );

    // Bls12-377 tests
    batch_field_test!(
        [
            test_lookup_fuctionality_1,
            test_lookup_fuctionality_2
        ],
        [] => bls12_377_scalar_field
    );
}
