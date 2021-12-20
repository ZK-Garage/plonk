// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-INFRA. All rights reserved.

//! Module containing the lookup functionality

use crate::error::Error;
use crate::lookup::MultiSet;
use ark_ff::Field;

/// This struct is a table, contaning a vector,
/// of arity 4 where each of the values is a
/// scalar. The elements of the table are
/// determined by the function g for
/// g(x,y), used to compute tuples.
///
/// This struct will be used to determine
/// the outputs of gates within arithmetic
/// circuits.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct LookupTable<F>(pub Vec<[F; 4]>)
where
    F: Field;

impl<F> LookupTable<F>
where
    F: Field,
{
    /// Create a new, empty Plookup table, with arity 4.
    pub fn new() -> Self {
        Default::default()
    }

    /// Insert a new row for an addition operation.
    /// This function needs to know the upper bound of the amount of addition
    /// operations that will be done in the plookup table.
    pub fn insert_add_row(&mut self, a: u64, b: u64, upper_bound: u64) {
        let c = (a + b) % upper_bound;
        self.0.push([F::from(a), F::from(b), F::from(c), F::zero()]);
    }

    /// Insert a new row for an addition operation.
    /// This function needs to know the upper bound of the amount of addition
    /// operations that will be done in the plookup table.
    pub fn insert_special_row(&mut self, a: F, b: F, c: F, d: F) {
        self.0.push([a, b, c, d]);
    }

    /// Insert a new row for an multiplication operation.
    /// This function needs to know the upper bound of the amount of
    /// multiplication operations that will be done in the plookup table.
    pub fn insert_mul_row(&mut self, a: u64, b: u64, upper_bound: u64) {
        let c = (a * b) % upper_bound;
        self.0.push([F::from(a), F::from(b), F::from(c), F::one()]);
    }

    /// Insert a new row for an XOR operation.
    /// This function needs to know the upper bound of the amount of XOR
    /// operations that will be done in the plookup table.
    pub fn insert_xor_row(&mut self, a: u64, b: u64, upper_bound: u64) {
        let c = (a ^ b) % upper_bound;
        self.0.push([F::from(a), F::from(b), F::from(c), -F::one()]);
    }

    /// Insert a new row for an AND operation.
    /// This function needs to know the upper bound of the amount of XOR
    /// operations that will be done in the plookup table.
    pub fn insert_and_row(&mut self, a: u64, b: u64, upper_bound: u64) {
        let c = (a & b) % upper_bound;
        self.0
            .push([F::from(a), F::from(b), F::from(c), F::from(2u64)]);
    }

    /// Function builds a table from more than one operation. This is denoted
    /// as 'Multiple Tables' in the paper. If, for example, we are using lookup
    /// tables for both XOR and mul operataions, we can create a table where the
    /// rows 0..n/2 use a mul function on all 2^n indices and have the 4th wire
    /// storing index 0. For all indices n/2..n, an XOR gate can be added, where
    /// the index of the 4th wire is 0.
    /// These numbers require exponentiation outside, for the lower bound,
    /// otherwise the range cannot start from zero, as 2^0 = 1.
    pub fn insert_multi_add(&mut self, lower_bound: u64, n: u8) {
        let upper_bound = 2u64.pow(n.into());
        let range = lower_bound..upper_bound;
        for a in range.clone() {
            range
                .clone()
                .for_each(|b| self.insert_add_row(a, b, upper_bound));
        }
    }

    /// Function builds a table from mutiple operations. If, for example,
    /// we are using lookup tables for both XOR and mul operataions, we can
    /// create a table where the rows 0..n/2 use a mul function on all 2^n
    /// indices and have the 4th wire storing index 0. For all indices n/2..n,
    /// an XOR gate can be added, wheren the index of the 4th wire is 0.
    /// These numbers require exponentiation outside, for the lower bound,
    /// otherwise the range cannot start from zero, as 2^0 = 1.
    /// Particular multiplication row(s) can be added with this function.
    pub fn insert_multi_mul(&mut self, lower_bound: u64, n: u8) {
        let upper_bound = 2u64.pow(n.into());
        let range = lower_bound..upper_bound;
        for a in range.clone() {
            range
                .clone()
                .for_each(|b| self.insert_mul_row(a, b, upper_bound));
        }
    }

    /// Function builds a table from mutiple operations. If, for example,
    /// we are using lookup tables for both XOR and mul operataions, we can
    /// create a table where the rows 0..n/2 use a mul function on all 2^n
    /// indices and have the 4th wire storing index 0. For all indices n/2..n,
    /// an XOR gate can be added, wheren the index of the 4th wire is 0.
    /// These numbers require exponentiation outside, for the lower bound,
    /// otherwise the range cannot start from zero, as 2^0 = 1.
    /// Particular XOR row(s) can be added with this function.
    pub fn insert_multi_xor(&mut self, lower_bound: u64, n: u8) {
        let upper_bound = 2u64.pow(n.into());
        let range = lower_bound..upper_bound;
        for a in range.clone() {
            range
                .clone()
                .for_each(|b| self.insert_xor_row(a, b, upper_bound));
        }
    }

    /// Function builds a table from mutiple operations. If, for example,
    /// we are using lookup tables for both XOR and mul operataions, we can
    /// create a table where the rows 0..n/2 use a mul function on all 2^n
    /// indices and have the 4th wire storing index 0. For all indices n/2..n,
    /// an XOR gate can be added, wheren the index of the 4th wire is 0.
    /// These numbers require exponentiation outside, for the lower bound,
    /// otherwise the range cannot start from zero, as 2^0 = 1.
    /// Particular AND row(s) can be added with this function.
    pub fn insert_multi_and(&mut self, lower_bound: u64, n: u8) {
        let upper_bound = 2u64.pow(n.into());
        let range = lower_bound..upper_bound;
        for a in range.clone() {
            range
                .clone()
                .for_each(|b| self.insert_and_row(a, b, upper_bound));
        }
    }

    /// Takes in a table, which is a list of vectors containing
    /// 4 elements, and turns them into 4 distinct multisets for
    /// a, b, c and d.
    pub fn vec_to_multiset(
        &self,
    ) -> (MultiSet<F>, MultiSet<F>, MultiSet<F>, MultiSet<F>) {
        let mut multiset_a = MultiSet::new();
        let mut multiset_b = MultiSet::new();
        let mut multiset_c = MultiSet::new();
        let mut multiset_d = MultiSet::new();
        self.0.iter().for_each(|row| {
            multiset_a.push(row[0]);
            multiset_b.push(row[1]);
            multiset_c.push(row[2]);
            multiset_d.push(row[3]);
        });
        (multiset_a, multiset_b, multiset_c, multiset_d)
    }

    /// Attempts to find an output value, given two input values, by querying
    /// the lookup table. The final wire holds the index of the table. The
    /// element must be predetermined to be between -1 and 2 depending on
    /// the type of table used. If the element does not exist, it will
    /// return an error.
    pub fn lookup(&self, a: F, b: F, d: F) -> Result<F, Error> {
        let pos = self
            .0
            .iter()
            .position(|row| row[0] == a && row[1] == b && row[3] == d)
            .ok_or(Error::ElementNotIndexed)?;

        Ok(self.0[pos][2])
    }
}

#[cfg(test)]
mod test {
    // TODO: use super::*;

    // FIXME: Run tests on both BLS fields.

    /* FIXME: Implement ADD table.
    fn test_add_table<F>()
    where
        F: Field,
    {
        let n = 4;

        let table = LookupTable::add_table(0, n);

        // Create an identical matrix, but with std numbers.
        // This way, we can also do the modulo operation, and properly
        // check all results.
        let mut i = 0;
        let p = 2u64.pow(n as u32);
        (0..p).for_each(|a| {
            (0..p).for_each(|b| {
                let c = (a + b) % p;
                assert_eq!(F::from(c), table.0[i][2]);
                i += 1;
            })
        });

        assert_eq!(
            table.0.len() as u64,
            2u64.pow(n as u32) * 2u64.pow(n as u32)
        );
    }
    */

    /* FIXME: Implement XOR table.
    fn test_xor_table<F>()
    where
        F: Field,
    {
        let n = 4;

        let table = LookupTable::xor_table(0, n);

        // println!("{:?}", table);
        let mut i = 0;
        let p = 2u64.pow(n as u32);
        (0..p).for_each(|a| {
            (0..p).for_each(|b| {
                let c = (a ^ b) % p;
                assert_eq!(F::from(c), table.0[i][2]);
                i += 1;
            })
        });

        assert_eq!(
            table.0.len() as u64,
            2u64.pow(n as u32) * 2u64.pow(n as u32)
        );
    }
    */

    /* FIXME: Implement MUL table.
    fn test_mul_table<F>()
    where
        F: Field,
    {
        let n = 4;

        let table = LookupTable::mul_table(0, n);

        // println!("{:?}", table);
        let mut i = 0;
        let p = 2u64.pow(n as u32);
        (0..p).for_each(|a| {
            (0..p).for_each(|b| {
                let c = (a * b) % p;
                assert_eq!(F::from(c), table.0[i][2]);
                i += 1;
            })
        });

        assert_eq!(
            table.0.len() as u64,
            2u64.pow(n as u32) * 2u64.pow(n as u32)
        );
    }
    */

    /* FIXME: Implement ADD table.
    fn test_lookup_arity_3<F>()
    where
        F: Field,
    {
        let add_table = LookupTable::add_table(0, 3);

        assert!(add_table.lookup(F::from(2u32), F::from(3u32)).is_ok());

        let output = add_table.0[1][0] + add_table.0[1][1];

        assert_eq!(output, F::one());

        let second_output = add_table.0[12][0] + add_table.0[12][1];

        assert_eq!(second_output, F::from(5u32));
    }
    */

    /* FIXME: Implement XOR table.
    fn test_missing_lookup_value<F>()
    where
        F: Field,
    {
        let xor_table = LookupTable::xor_table(0, 5);
        assert!(xor_table.lookup(F::from(17u32), F::from(367u32)).is_err());
    }
    */

    /* FIXME: Check XOR implemenation for fields.
    fn test_concatenated_table<F>()
    where
        F: Field,
    {
        let mut table = LookupTable::<F>::new();

        table.insert_multi_xor(0, 5);
        table.insert_multi_add(4, 7);

        assert_eq!(table.0.last().unwrap()[2], F::from(126u64));
        let xor: F = table.0[36][0] ^ table.0[36][1];
        assert_eq!(xor, F::from(5u64));
    }
    */
}