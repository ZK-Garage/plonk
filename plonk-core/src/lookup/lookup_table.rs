// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::error::Error;
use crate::lookup::MultiSet;
use ark_ff::Field;

/// This struct is a table, contaning a vector, of arity 4 where each of the
/// values is a scalar. The elements of the table are determined by the function
/// g for g(x,y), used to compute tuples.
///
/// This struct will be used to determine the outputs of gates within arithmetic
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

    /// Returns the length of the `LookupTable` vector.
    pub fn size(&self) -> usize {
        self.0.len()
    }

    /// Pushes a row to the `LookupTable` vector.
    fn push(&mut self, row: [F; 4]) {
        self.0.push(row);
    }

    /// Insert a new row
    pub fn insert_row(&mut self, a: F, b: F, c: F, d: F) {
        self.push([a, b, c, d]);
    }

    /// Insert a new row for an addition operation.
    ///
    /// This function needs to know the upper bound of the amount of addition
    /// operations that will be done in the plookup table.
    /// The result will be:  a + b mod 2^uppder_bound
    pub fn insert_add_row(&mut self, a: u64, b: u64, upper_bound: u64) {
        let c = (a + b) % upper_bound;
        self.insert_row(F::from(a), F::from(b), F::from(c), F::zero());
    }

    /// Insert a new row for a multiplication operation.
    ///
    /// This function needs to know the upper bound of the amount of
    /// multiplication operations that will be done in the plookup table.
    /// The result will be:  a * b mod 2^uppder_bound
    pub fn insert_mul_row(&mut self, a: u64, b: u64, upper_bound: u64) {
        let c = (a * b) % upper_bound;
        self.insert_row(F::from(a), F::from(b), F::from(c), F::one());
    }

    /// Insert a new row for an XOR operation.
    ///
    /// This function needs to know the upper bound of the amount of XOR
    /// operations that will be done in the plookup table.
    /// The result will be:  a XOR b mod 2^uppder_bound
    pub fn insert_xor_row(&mut self, a: u64, b: u64, upper_bound: u64) {
        let c = (a ^ b) % upper_bound;
        self.insert_row(F::from(a), F::from(b), F::from(c), -F::one());
    }

    /// Insert a new row for an AND operation.
    ///
    /// This function needs to know the upper bound of the amount of AND
    /// operations that will be done in the plookup table.
    /// The result will be:  a AND b mod 2^uppder_bound
    pub fn insert_and_row(&mut self, a: u64, b: u64, upper_bound: u64) {
        let c = (a & b) % upper_bound;
        self.insert_row(F::from(a), F::from(b), F::from(c), F::from(2u64));
    }

    /// Function builds a table from more than one operation. This is denoted
    /// as 'Multiple Tables' in the paper. If, for example, we are using lookup
    /// tables for both XOR and mul operataions, we can create a table where the
    /// rows 0..n/2 use a mul function on all 2^n indices and have the 4th wire
    /// storing index 0. For all indices n/2..n, an XOR gate can be added, where
    /// the index of the 4th wire is 0.
    /// These numbers require exponentiation outside, for the lower bound,
    /// otherwise the range cannot start from zero, as 2^0 = 1.
    pub fn insert_multi_add(&mut self, lower_bound: u64, n: u32) {
        let upper_bound = 2u64.pow(n);
        for a in lower_bound..upper_bound {
            for b in lower_bound..upper_bound {
                self.insert_add_row(a, b, upper_bound);
            }
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
    pub fn insert_multi_mul(&mut self, lower_bound: u64, n: u32) {
        let upper_bound = 2u64.pow(n);
        for a in lower_bound..upper_bound {
            for b in lower_bound..upper_bound {
                self.insert_mul_row(a, b, upper_bound);
            }
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
    pub fn insert_multi_xor(&mut self, lower_bound: u64, n: u32) {
        let upper_bound = 2u64.pow(n);
        for a in lower_bound..upper_bound {
            for b in lower_bound..upper_bound {
                self.insert_xor_row(a, b, upper_bound);
            }
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
    pub fn insert_multi_and(&mut self, lower_bound: u64, n: u32) {
        let upper_bound = 2u64.pow(n);
        for a in lower_bound..upper_bound {
            for b in lower_bound..upper_bound {
                self.insert_and_row(a, b, upper_bound);
            }
        }
    }

    /// Takes in a table, which is a vector of slices containing
    /// 4 elements, and turns them into 4 distinct multisets for
    /// a, b, c and d.
    pub fn vec_to_multiset(&self) -> Vec<MultiSet<F>> {
        let mut result = vec![MultiSet::new(); 4];
        self.0.iter().for_each(|row| {
            result.iter_mut().enumerate().for_each(|(index, multiset)| {
                multiset.push(row[index]);
            })
        });
        result
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

    /// Creates an addition table for addends from the lower bound up to the
    /// upper bound 2^n
    pub fn add_table(lower_bound: u64, n: u32) -> LookupTable<F> {
        let mut table = LookupTable::new();
        table.insert_multi_add(lower_bound, n);
        table
    }

    /// Creates an xor table for addends from the lower bound up to the upper
    /// bound 2^n
    pub fn xor_table(lower_bound: u64, n: u32) -> LookupTable<F> {
        let mut table = LookupTable::new();
        table.insert_multi_xor(lower_bound, n);
        table
    }

    /// Creates an addition table for addends from the lower bound up to the
    /// upper bound 2^n
    pub fn mul_table(lower_bound: u64, n: u32) -> LookupTable<F> {
        let mut table = LookupTable::new();
        table.insert_multi_mul(lower_bound, n);
        table
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::batch_field_test;
    use ark_bls12_377::Fr as Bls12_377_scalar_field;
    use ark_bls12_381::Fr as Bls12_381_scalar_field;

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
        let p = 2u64.pow(n);
        (0..p).for_each(|a| {
            (0..p).for_each(|b| {
                let c = (a + b) % p;
                assert_eq!(F::from(c), table.0[i][2]);
                i += 1;
            })
        });
        assert_eq!(table.0.len() as u64, 2u64.pow(n) * 2u64.pow(n));
    }

    fn test_xor_table<F>()
    where
        F: Field,
    {
        let n = 4;
        let table = LookupTable::xor_table(0, n);
        let mut i = 0;
        let p = 2u64.pow(n);
        (0..p).for_each(|a| {
            (0..p).for_each(|b| {
                let c = a ^ b;
                assert_eq!(F::from(c), table.0[i][2]);
                i += 1;
            })
        });
        assert_eq!(table.0.len() as u64, 2u64.pow(n) * 2u64.pow(n));
    }

    fn test_mul_table<F>()
    where
        F: Field,
    {
        let n = 4;
        let table = LookupTable::mul_table(0, n);
        let mut i = 0;
        let p = 2u64.pow(n);
        (0..p).for_each(|a| {
            (0..p).for_each(|b| {
                let c = (a * b) % p;
                assert_eq!(F::from(c), table.0[i][2]);
                i += 1;
            })
        });
        assert_eq!(table.0.len() as u64, 2u64.pow(n) * 2u64.pow(n));
    }

    fn test_lookup_arity_3<F>()
    where
        F: Field,
    {
        let add_table = LookupTable::add_table(0, 3);
        assert!(add_table
            .lookup(F::from(2u32), F::from(3u32), F::from(0u32))
            .is_ok());
        let output = add_table.0[1][0] + add_table.0[1][1];
        assert_eq!(output, F::one());
        let second_output = add_table.0[12][0] + add_table.0[12][1];
        assert_eq!(second_output, F::from(5u32));
    }

    fn test_missing_lookup_value<F>()
    where
        F: Field,
    {
        let xor_table = LookupTable::xor_table(0, 5);
        assert!(xor_table
            .lookup(F::from(17u32), F::from(367u32), F::from(0u32))
            .is_err());
    }

    fn test_concatenated_table<F: ark_ff::PrimeField>()
    where
        F: Field,
    {
        use ark_ff::BigInteger;
        let mut table = LookupTable::<F>::new();
        table.insert_multi_xor(0, 5);
        table.insert_multi_add(4, 7);
        assert_eq!(table.0.last().unwrap()[2], F::from(126u64));
        let xor: F = F::from_repr(BigInteger::from_bits_le(
            table.0[36][0]
                .into_repr()
                .to_bits_le()
                .iter()
                .zip(table.0[36][1].into_repr().to_bits_le().iter())
                .map(|(l, r)| l ^ r)
                .collect::<Vec<_>>()
                .as_slice(),
        ))
        .unwrap();

        assert_eq!(xor, F::from(5u64));
    }

    // Bls12-381 tests
    batch_field_test!(
        [
            test_add_table,
            test_xor_table,
            test_mul_table,
            test_lookup_arity_3,
            test_missing_lookup_value,
            test_concatenated_table
        ],
        [] => Bls12_381_scalar_field
    );

    // Bls12-377 tests
    batch_field_test!(
        [
            test_add_table,
            test_xor_table,
            test_mul_table,
            test_lookup_arity_3,
            test_missing_lookup_value,
            test_concatenated_table
        ],
        [] => Bls12_377_scalar_field
    );
}
