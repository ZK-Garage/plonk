// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Boolean Gates

use crate::constraint_system::{StandardComposer, Variable};
use ark_ec::ModelParameters;
use ark_ff::PrimeField;

impl<F, P> StandardComposer<F, P>
where
    F: PrimeField,
    P: ModelParameters<BaseField = F>,
{
    /// Adds a boolean constraint (also known as binary constraint) where
    /// the gate eq. will enforce that the [`Variable`] received is either `0`
    /// or `1` by adding a constraint in the circuit.
    ///
    /// Note that using this constraint with whatever [`Variable`] that is not
    /// representing a value equalling 0 or 1, will always force the equation to
    /// fail.
    pub fn boolean_gate(&mut self, a: Variable) -> Variable {
        self.w_l.push(a);
        self.w_r.push(a);
        self.w_o.push(a);
        self.w_4.push(self.zero_var);

        self.q_m.push(F::one());
        self.q_l.push(F::zero());
        self.q_r.push(F::zero());
        self.q_o.push(-F::one());
        self.q_c.push(F::zero());
        self.q_4.push(F::zero());
        self.q_arith.push(F::one());

        // add high degree selectors
        self.q_hl.push(F::zero());
        self.q_hr.push(F::zero());
        self.q_h4.push(F::zero());

        self.q_range.push(F::zero());
        self.q_logic.push(F::zero());
        self.q_fixed_group_add.push(F::zero());
        self.q_variable_group_add.push(F::zero());
        self.q_lookup.push(F::zero());

        self.perm
            .add_variables_to_map(a, a, a, self.zero_var, self.n);

        self.n += 1;

        a
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
    use ark_ec::TEModelParameters;

    fn test_correct_bool_gate<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                let zero = composer.zero_var();
                let one = composer.add_input(F::one());
                composer.boolean_gate(zero);
                composer.boolean_gate(one);
            },
            32,
        );
        assert!(res.is_ok())
    }

    fn test_incorrect_bool_gate<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                let zero = composer.add_input(F::from(5u64));
                let one = composer.add_input(F::one());
                composer.boolean_gate(zero);
                composer.boolean_gate(one);
            },
            32,
        );
        assert!(res.is_err())
    }

    // Test for Bls12_381
    batch_test!(
        [
            test_correct_bool_gate,
            test_incorrect_bool_gate
        ],
        [] => (
            Bls12_381, ark_ed_on_bls12_381::EdwardsParameters
        )
    );

    // Test for Bls12_377
    batch_test!(
        [
            test_correct_bool_gate,
            test_incorrect_bool_gate
        ],
        [] => (
            Bls12_377, ark_ed_on_bls12_377::EdwardsParameters        )
    );
}
