// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-INFRA. All rights reserved.

use crate::error::Error;
use crate::lookup::{LookupTable, MultiSet};
use ark_ff::Field;

/// This witness table contains quieries
/// to a lookup table for lookup gates
/// This table is of arity 3.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct WitnessTable<F>
where
    F: Field,
{
    /// This column represents the
    /// first values inside the lookup
    /// table. At gate checks, this
    /// can be regarded as the first
    /// wire
    pub f_1: MultiSet<F>,

    /// This column represents the
    /// first values inside the lookup
    /// table. At gate checks, this
    /// can be regarded as the second
    /// wire
    pub f_2: MultiSet<F>,

    /// This column represents the
    /// first values inside the lookup
    /// table. At gate checks, this
    /// can be regarded as the third
    /// wire
    pub f_3: MultiSet<F>,

    /// This column represents the
    /// first values inside the lookup
    /// table. At gate checks, this
    /// can be regarded as the fourth
    /// wire
    pub f_4: MultiSet<F>,
}

impl<F> WitnessTable<F>
where
    F: Field,
{
    /// Initialses empty witness table of arity 4
    pub fn new() -> Self {
        Default::default()
    }

    /// This allows the witness table to be filled directly without
    /// taking any vaules, or the the results, from the lookup table.
    /// If the values do no exists in the lookup table, then the proof
    /// will fail when witness and preprocessed tables are concatenated.
    pub fn from_wire_values(
        &mut self,
        left_wire_val: F,
        right_wire_val: F,
        output_wire_val: F,
        fourth_wire_val: F,
    ) {
        self.f_1.push(left_wire_val);
        self.f_2.push(right_wire_val);
        self.f_3.push(output_wire_val);
        self.f_4.push(fourth_wire_val);
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
        let output_wire_val = lookup_table.lookup(
            left_wire_val,
            right_wire_val,
            fourth_wire_val,
        )?;
        self.f_1.push(left_wire_val);
        self.f_2.push(right_wire_val);
        self.f_3.push(output_wire_val);
        self.f_4.push(fourth_wire_val);
        Ok(())
    }
}
#[cfg(test)]
mod test {
    use super::*;
    use crate::lookup::LookupTable;

    // FIXME: Run tests on both BLS fields.

    /* FIXME: Implement XOR table.
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
            .value_from_table(&lookup_table, F::from(2), F::from(5))
            .is_ok());
        // Check that non existent elements cause a failure
        assert!(f
            .value_from_table(&lookup_table, F::from(25), F::from(5))
            .is_err());
    }
    */

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
}