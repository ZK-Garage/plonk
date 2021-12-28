// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-INFRA. All rights reserved.

//! Boolean Gates

use crate::constraint_system::{StandardComposer, Variable};
use ark_ec::{PairingEngine, TEModelParameters};
use ark_ff::Field;
use num_traits::{One, Zero};

use super::arithmetic::ArithmeticGate;

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

    /// A gate which outputs a variable equal to 0 if the value
    /// of its input is zero and outputs 1 otherwise.
    pub fn is_nonzero_gate(&mut self, a: Variable) -> Variable {
        //Get inverse value or just give 1 here if a is zero
        let a_inv_value = self.variables.get(&a).unwrap().inverse().unwrap_or(E::Fr::one());
        let a_inv = self.add_input(a_inv_value);
        
        // This variable has value zero if input a is zero and value 1 otherwise
        let product = self.arithmetic_gate(|gate| gate.witness(a, a_inv, None).mul(E::Fr::one()));
        product
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

    fn test_correct_is_nonzero_gate<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let one = composer.add_input(E::Fr::one());
                let two = composer.add_input(E::Fr::one().double());
                let is_nonzero = composer.is_nonzero_gate(two);
                composer.assert_equal(is_nonzero, one);
            },
             32
            );
            assert!(res.is_ok())
    }

    fn test_correct_is_nonzero_gate2<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let zero = composer.zero_var();
                let is_nonzero = composer.is_nonzero_gate(zero);
                composer.assert_equal(is_nonzero, zero);
            },
             32
            );
            assert!(res.is_ok())
    }

    fn test_incorrect_is_nonzero_gate<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let two = composer.add_input(E::Fr::one().double());
                let is_nonzero = composer.is_nonzero_gate(two);
                composer.assert_equal(is_nonzero, composer.zero_var());
            },
             32
            );
            assert!(res.is_err())
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
            test_correct_is_nonzero_gate,
            test_correct_is_nonzero_gate2,
            test_incorrect_is_nonzero_gate,
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
            test_correct_is_nonzero_gate,
            test_correct_is_nonzero_gate2,
            test_incorrect_is_nonzero_gate,
            test_correct_bool_gate,
            test_incorrect_bool_gate
        ],
        [] => (
            Bls12_377,
            ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
