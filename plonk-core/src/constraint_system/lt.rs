// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) ZK-GARAGE . All rights reserved.

//! Comparator gadget

use crate::constraint_system::{StandardComposer, Variable};
use ark_ec::TEModelParameters;
use ark_ff::PrimeField;

impl<F, P> StandardComposer<F, P>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    /// Compares 2 variables withing a range, returns `1` if lhs < rhs, 0
    /// otherwise.
    ///
    /// Both `rhs` and `lhs` must be `<= 2**n_bits`
    /// The number of bits `n_bits` must be even.
    pub fn less_than(
        &mut self,
        lhs: Variable,
        rhs: Variable,
        n_bits: usize,
    ) -> Variable {
        let lhs_val = self.variables[&lhs];
        let rhs_val = self.variables[&rhs];

        let range_val = F::from(2u64.pow(n_bits as u32));

        let (lt, diff) = if lhs_val < rhs_val {
            let lt = self.add_input(F::one());
            let diff = self.add_input(range_val + lhs_val - rhs_val);
            (lt, diff)
        } else {
            let lt = self.add_input(F::zero());
            let diff = self.add_input(lhs_val - rhs_val);
            (lt, diff)
        };

        // Constraints
        // 1. lt must be binary
        self.boolean_gate(lt);

        //2. lhs - rhs = diff - (lt * range)
        //  lhs + (-1) rhs + (-1) diff + (range) lt = 0
        self.arithmetic_gate(|gate| {
            gate.witness(lhs, rhs, Some(lt))
                .add(F::one(), -F::one())
                .fan_in_3(-F::one(), diff)
                .out(range_val)
        });

        // 3. Range check lhs, rhs and diff
        self.range_gate(lhs, n_bits);
        self.range_gate(rhs, n_bits);
        self.range_gate(diff, n_bits);

        lt
    }
}
#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        batch_test, commitment::HomomorphicCommitment,
        constraint_system::helper::*,
    };
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;

    fn test_less_than<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                let four = composer.add_input(F::from(4u64));
                let five = composer.add_input(F::from(5u64));

                let output = composer.less_than(four, five, 4);
                composer.constrain_to_constant(output, F::one(), None);

                let output = composer.less_than(five, four, 4);
                composer.constrain_to_constant(output, F::zero(), None);

                let output = composer.less_than(five, five, 4);
                composer.constrain_to_constant(output, F::zero(), None);
            },
            200,
        );
        assert!(res.is_ok(), "{:?}", res.err().unwrap());
    }

    batch_test!(
        [
            test_less_than
        ],
        [] => (
            Bls12_377, ark_ed_on_bls12_377::EdwardsParameters
        )
    );

    batch_test!(
        [
            test_less_than
        ],
        [] => (
            Bls12_381, ark_ed_on_bls12_381::EdwardsParameters
        )
    );
}
