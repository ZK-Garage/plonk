// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Boolean Gates

use crate::{
    constraint_system::{StandardComposer, Variable},
    parameters::CircuitParameters,
};
use ark_ff::{One, Zero};

impl<P> StandardComposer<P>
where
    P: CircuitParameters,
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

        self.q_m.push(P::ScalarField::one());
        self.q_l.push(P::ScalarField::zero());
        self.q_r.push(P::ScalarField::zero());
        self.q_o.push(-P::ScalarField::one());
        self.q_c.push(P::ScalarField::zero());
        self.q_4.push(P::ScalarField::zero());
        self.q_arith.push(P::ScalarField::one());

        self.q_range.push(P::ScalarField::zero());
        self.q_logic.push(P::ScalarField::zero());
        self.q_fixed_group_add.push(P::ScalarField::zero());
        self.q_variable_group_add.push(P::ScalarField::zero());
        self.q_lookup.push(P::ScalarField::zero());

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
        batch_test, constraint_system::helper::*, parameters::test::*,
    };

    fn test_correct_bool_gate<P>()
    where
        P: CircuitParameters,
    {
        let res = gadget_tester::<P>(
            |composer: &mut StandardComposer<P>| {
                let zero = composer.zero_var();
                let one = composer.add_input(P::ScalarField::one());
                composer.boolean_gate(zero);
                composer.boolean_gate(one);
            },
            32,
        );
        assert!(res.is_ok())
    }

    fn test_incorrect_bool_gate<P>()
    where
        P: CircuitParameters,
    {
        let res = gadget_tester::<P>(
            |composer: &mut StandardComposer<P>| {
                let zero = composer.add_input(P::ScalarField::from(5u64));
                let one = composer.add_input(P::ScalarField::one());
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
        [] => [
            Bls12_381_KZG, Bls12_381_IPA, Bls12_377_KZG, Bls12_377_IPA
        ]
    );
}
