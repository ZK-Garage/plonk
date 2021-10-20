// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::constraint_system::StandardComposer;
use crate::constraint_system::Variable;
use ark_ec::{PairingEngine, ProjectiveCurve, TEModelParameters};
use num_traits::{One, Zero};

impl<
        E: PairingEngine,
        T: ProjectiveCurve<BaseField = E::Fr>,
        P: TEModelParameters<BaseField = E::Fr>,
    > StandardComposer<E, T, P>
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

        self.q_m.push(E::Fr::one());
        self.q_l.push(E::Fr::zero());
        self.q_r.push(E::Fr::zero());
        self.q_o.push(-E::Fr::one());
        self.q_c.push(E::Fr::zero());
        self.q_4.push(E::Fr::zero());
        self.q_arith.push(E::Fr::one());

        self.q_range.push(E::Fr::zero());
        self.q_logic.push(E::Fr::zero());
        self.q_fixed_group_add.push(E::Fr::zero());
        self.q_variable_group_add.push(E::Fr::zero());

        self.perm
            .add_variables_to_map(a, a, a, self.zero_var, self.n);

        self.n += 1;

        a
    }
}

#[cfg(test)]
mod boolean_gates_tests {
    use crate::constraint_system::helper::*;
    use crate::constraint_system::StandardComposer;

    use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
    use ark_ed_on_bls12_381::{
        EdwardsParameters as JubjubParameters,
        EdwardsProjective as JubjubProjective,
    };
    use num_traits::One;
    #[test]
    fn test_correct_bool_gate() {
        let res = gadget_tester(
            |composer: &mut StandardComposer<
                Bls12_381,
                JubjubProjective,
                JubjubParameters,
            >| {
                let zero = composer.zero_var();
                let one = composer.add_input(BlsScalar::one());

                composer.boolean_gate(zero);
                composer.boolean_gate(one);
            },
            32,
        );
        assert!(res.is_ok())
    }

    #[test]
    fn test_incorrect_bool_gate() {
        let res = gadget_tester(
            |composer: &mut StandardComposer<
                Bls12_381,
                JubjubProjective,
                JubjubParameters,
            >| {
                let zero = composer.add_input(BlsScalar::from(5));
                let one = composer.add_input(BlsScalar::one());

                composer.boolean_gate(zero);
                composer.boolean_gate(one);
            },
            32,
        );
        assert!(res.is_err())
    }
}
