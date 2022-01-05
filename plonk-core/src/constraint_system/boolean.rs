// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Boolean Gates

use crate::constraint_system::{StandardComposer, Variable};
use ark_ec::{PairingEngine, TEModelParameters};
use ark_ff::Field;
use num_traits::{One, Zero};

impl<E, P> StandardComposer<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
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

    /// A gate which outputs a variable whose value is 1 if
    /// the input is 0 and whose value is 0 otherwise
    pub fn is_zero_gate(&mut self, a: Variable) -> Variable {
        //Get relevant field values
        let a_value = self.variables.get(&a).unwrap();
        let y_value = a_value.inverse().unwrap_or_else(E::Fr::one);
        //This has value 1 if input value is zero, value 0 otherwise
        let b_value = E::Fr::one() - *a_value * y_value;

        let y = self.add_input(y_value);
        let b = self.add_input(b_value);
        let zero = self.zero_var();

        // Enforce constraints. The constraint system being used here is
        // a * y + b - 1 = 0
        // a * b = 0
        // where y is auxiliary and b is the boolean (a == 0).
        let _a_times_b = self.arithmetic_gate(|gate| {
            gate.witness(a, b, Some(zero)).mul(E::Fr::one())
        });

        let _constraint1 = self.arithmetic_gate(|gate| {
            gate.witness(a, y, Some(zero))
                .mul(E::Fr::one())
                .fan_in_3(E::Fr::one(), b)
                .constant(-E::Fr::one())
        });

        b
    }

    /// A gate which outputs a variable whose value is 1 if the
    /// two input variables have equal values and whose value is 0 otherwise.
    pub fn is_eq_gate(&mut self, a: Variable, b: Variable) -> Variable {
        let difference = self.arithmetic_gate(|gate| {
            gate.witness(a, b, None).add(E::Fr::one(), -E::Fr::one())
        });
        // TODO: These dummy constraints are a hack to avoid a 2^k issue,
        // specifically the q_arith poly being constant.
        // It should be safe to remove this dummy constraint once that issue
        // is resolved.
        self.add_dummy_constraints();
        self.is_zero_gate(difference)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::batch_test;
    use crate::constraint_system::helper::*;
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use num_traits::One;

    fn test_correct_is_zero_gate<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        // Check that it gives true on zero input:
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let one = composer.add_input(E::Fr::one());
                let is_zero = composer.is_zero_gate(composer.zero_var());
                composer.assert_equal(is_zero, one);
            },
            32,
        );

        // Check that it gives false on non-zero input:
        let res2 = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let one = composer.add_input(E::Fr::one());
                let is_zero = composer.is_zero_gate(one);
                composer.assert_equal(is_zero, composer.zero_var());
            },
            32,
        );

        assert!(res.is_ok() && res2.is_ok())
    }

    fn test_correct_is_eq_gate<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        // Check that it gives true on equal inputs:
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let one = composer.add_input(E::Fr::one());

                let field_element = E::Fr::one().double();
                let a = composer.add_input(field_element);
                let b = composer.add_input(field_element);
                let is_eq = composer.is_eq_gate(a, b);
                composer.assert_equal(is_eq, one);
            },
            32,
        );

        // Check that it gives false on non-equal inputs:
        let res2 = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let field_element = E::Fr::one().double();
                let a = composer.add_input(field_element);
                let b = composer.add_input(field_element.double());
                let is_eq = composer.is_eq_gate(a, b);
                composer.assert_equal(is_eq, composer.zero_var());
            },
            32,
        );

        assert!(res.is_ok() && res2.is_ok())
    }

    fn test_correct_bool_gate<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let zero = composer.zero_var();
                let one = composer.add_input(E::Fr::one());
                composer.boolean_gate(zero);
                composer.boolean_gate(one);
            },
            32,
        );
        assert!(res.is_ok())
    }

    fn test_incorrect_bool_gate<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let zero = composer.add_input(E::Fr::from(5u64));
                let one = composer.add_input(E::Fr::one());
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
            test_correct_is_zero_gate,
            test_correct_is_eq_gate,
            test_correct_bool_gate,
            test_incorrect_bool_gate
        ],
        [] => (
            Bls12_381,
            ark_ed_on_bls12_381::EdwardsParameters
        )
    );

    // Test for Bls12_377
    batch_test!(
        [
            test_correct_is_zero_gate,
            test_correct_is_eq_gate,
            test_correct_bool_gate,
            test_incorrect_bool_gate
        ],
        [] => (
            Bls12_377,
            ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
