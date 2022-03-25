// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) ZK-Garage. All rights reserved.

use crate::constraint_system::{StandardComposer, Variable};
use ark_ec::TEModelParameters;
use ark_ff::PrimeField;

impl<F, P> StandardComposer<F, P>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    /// Adds a plookup gate to the circuit with its corresponding
    /// constraints.
    ///
    /// This type of gate is usually used when we need to have
    /// the largest amount of performance and the minimum circuit-size
    /// possible. Since it allows the end-user to set every selector coefficient
    /// as scaling value on the gate eq.
    pub fn lookup_gate(
        &mut self,
        a: Variable,
        b: Variable,
        c: Variable,
        d: Option<Variable>,
        q_l: F,
        q_r: F,
        q_o: F,
        q_4: F,
        q_c: F,
        pi: Option<F>,
    ) -> Variable {
        // Check if advice wire has a value
        let d = match d {
            Some(var) => var,
            None => self.zero_var,
        };

        self.w_l.push(a);
        self.w_r.push(b);
        self.w_o.push(c);
        self.w_4.push(d);

        // Add selector vectors
        self.q_l.push(q_l);
        self.q_r.push(q_r);
        self.q_o.push(q_o);
        self.q_c.push(q_c);
        self.q_4.push(q_4);
        self.q_arith.push(F::zero());
        self.q_m.push(F::zero());
        self.q_range.push(F::zero());
        self.q_logic.push(F::zero());
        self.q_fixed_group_add.push(F::zero());
        self.q_variable_group_add.push(F::zero());

        // For a lookup gate, only one selector poly is
        // turned on as the output is inputted directly
        self.q_lookup.push(F::one());

        if let Some(pi) = pi {
            assert!(self
                .public_inputs_sparse_store
                .insert(self.n, pi)
                .is_none());
        }

        self.perm.add_variables_to_map(a, b, c, d, self.n);

        self.n += 1;

        c
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        batch_test, commitment::HomomorphicCommitment,
        constraint_system::helper::*, lookup::LookupTable,
    };
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;

    fn test_plookup_xor<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                let xor_1bit = LookupTable::<F>::xor_table(0, 1);
                composer.lookup_table = xor_1bit;
                let zero = composer.zero_var;
                let one = composer.add_input(F::one());
                let negative_one = composer.add_input(-F::one());

                composer.lookup_gate(
                    one,
                    one,
                    zero,
                    Some(negative_one),
                    F::zero(),
                    F::zero(),
                    F::zero(),
                    F::zero(),
                    F::zero(),
                    None,
                );
            },
            4,
        );
        assert!(res.is_ok(), "{:?}", res.err().unwrap());
    }

    // Bls12-381 tests
    batch_test!(
        [
            test_plookup_xor
        ],
        [] => (
            Bls12_381, ark_ed_on_bls12_381::EdwardsParameters
        )
    );

    // Bls12-377 tests
    batch_test!(
        [
            test_plookup_xor
        ],
        [] => (
            Bls12_377, ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
