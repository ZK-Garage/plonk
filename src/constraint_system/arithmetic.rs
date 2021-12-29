// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Simple Arithmetic Gates

use crate::constraint_system::StandardComposer;
use crate::constraint_system::Variable;
use ark_ec::{PairingEngine, TEModelParameters};
use num_traits::{One, Zero};

#[derive(Debug, Clone, Copy)]
pub struct ArithmeticGate<E: PairingEngine> {
    pub(crate) witness: Option<(Variable, Variable, Option<Variable>)>,
    pub(crate) fan_in_3: Option<(E::Fr, Variable)>,
    pub(crate) mul_selector: E::Fr,
    pub(crate) add_selectors: (E::Fr, E::Fr),
    pub(crate) out_selector: E::Fr,
    pub(crate) const_selector: E::Fr,
    pub(crate) pi: Option<E::Fr>,
}

impl<E: PairingEngine> Default for ArithmeticGate<E> {
    fn default() -> Self {
        Self {
            witness: None,
            fan_in_3: None,
            mul_selector: E::Fr::zero(),
            add_selectors: (E::Fr::zero(), E::Fr::zero()),
            out_selector: -E::Fr::one(),
            const_selector: E::Fr::zero(),
            pi: None,
        }
    }
}

impl<E: PairingEngine> ArithmeticGate<E> {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn witness(
        &mut self,
        w_l: Variable,
        w_r: Variable,
        w_o: Option<Variable>,
    ) -> &mut Self {
        self.witness = Some((w_l, w_r, w_o));
        self
    }

    pub fn fan_in_3(&mut self, q_4: E::Fr, w_4: Variable) -> &mut Self {
        self.fan_in_3 = Some((q_4, w_4));
        self
    }

    pub fn mul(&mut self, q_m: E::Fr) -> &mut Self {
        self.mul_selector = q_m;
        self
    }

    pub fn add(&mut self, q_l: E::Fr, q_r: E::Fr) -> &mut Self {
        self.add_selectors = (q_l, q_r);
        self
    }

    pub fn out(&mut self, q_o: E::Fr) -> &mut Self {
        self.out_selector = q_o;
        self
    }

    pub fn constant(&mut self, q_c: E::Fr) -> &mut Self {
        self.const_selector = q_c;
        self
    }

    pub fn pi(&mut self, pi: E::Fr) -> &mut Self {
        self.pi = Some(pi);
        self
    }

    pub fn build(&mut self) -> Self {
        *self
    }
}

impl<E: PairingEngine, P: TEModelParameters<BaseField = E::Fr>>
    StandardComposer<E, P>
{
    /// Function used to generate any arithmetic gate with fan-in-2 or fan-in-3.
    pub fn arithmetic_gate<F>(&mut self, func: F) -> Variable
    where
        F: FnOnce(&mut ArithmeticGate<E>) -> &mut ArithmeticGate<E>,
    {
        let gate = {
            let mut gate = ArithmeticGate::<E>::new();
            func(&mut gate).build()
        };

        if gate.witness.is_none() {
            panic!("Missing left and right wire witnesses")
        }

        let (w4, q4) = if gate.fan_in_3.is_none() {
            self.w_4.push(self.zero_var);
            self.q_4.push(E::Fr::zero());
            (self.zero_var, E::Fr::zero())
        } else {
            let (q4, w4) = gate.fan_in_3.unwrap();
            self.w_4.push(w4);
            self.q_4.push(q4);
            (w4, q4)
        };

        let gate_witness = gate.witness.unwrap();
        self.w_l.push(gate_witness.0);
        self.w_r.push(gate_witness.1);
        self.q_l.push(gate.add_selectors.0);
        self.q_r.push(gate.add_selectors.1);

        // Add selector vectors
        self.q_m.push(gate.mul_selector);
        self.q_o.push(gate.out_selector);
        self.q_c.push(gate.const_selector);

        self.q_arith.push(E::Fr::one());
        self.q_range.push(E::Fr::zero());
        self.q_logic.push(E::Fr::zero());
        self.q_fixed_group_add.push(E::Fr::zero());
        self.q_variable_group_add.push(E::Fr::zero());

        if let Some(pi) = gate.pi {
            let insert_res = self.public_inputs_sparse_store.insert(self.n, pi);
            assert!(
                insert_res.is_none(),
                "Attempting to overwrite an already existing PI"
            )
        };

        let c = gate_witness.2.unwrap_or_else(|| {
            self.add_input(
                ((gate.mul_selector
                    * (self.variables[&gate_witness.0]
                        * self.variables[&gate_witness.1]))
                    + gate.add_selectors.0 * self.variables[&gate_witness.0]
                    + gate.add_selectors.1 * self.variables[&gate_witness.1]
                    + gate.const_selector
                    + q4 * self.variables[&w4]
                    + gate.pi.unwrap_or_default())
                    * (-gate.out_selector),
            )
        });
        self.w_o.push(c);
        self.perm.add_variables_to_map(
            gate_witness.0,
            gate_witness.1,
            c,
            w4,
            self.n,
        );
        self.n += 1;

        c
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::batch_test;
    use crate::constraint_system::helper::*;
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_ec::PairingEngine;
    use ark_ec::TEModelParameters;
    use num_traits::One;

    fn test_public_inputs<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let var_one = composer.add_input(E::Fr::one());

                let should_be_three = composer.arithmetic_gate(|gate| {
                    gate.witness(var_one, var_one, None)
                        .add(E::Fr::one(), E::Fr::one())
                        .pi(E::Fr::one())
                });

                composer.constrain_to_constant(
                    should_be_three,
                    E::Fr::from(3u64),
                    None,
                );

                let should_be_four = composer.arithmetic_gate(|gate| {
                    gate.witness(var_one, var_one, None)
                        .add(E::Fr::one(), E::Fr::one())
                        .pi(E::Fr::from(2u64))
                });

                composer.constrain_to_constant(
                    should_be_four,
                    E::Fr::from(4u64),
                    None,
                );
            },
            200,
        );
        assert!(res.is_ok(), "{:?}", res.err().unwrap());
    }

    fn test_correct_add_mul_gate<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                // Verify that (4+5+5) * (6+7+7) = 280
                let four = composer.add_input(E::Fr::from(4u64));
                let five = composer.add_input(E::Fr::from(5u64));
                let six = composer.add_input(E::Fr::from(6u64));
                let seven = composer.add_input(E::Fr::from(7u64));

                let fourteen = composer.arithmetic_gate(|gate| {
                    gate.witness(four, five, None)
                        .add(E::Fr::one(), E::Fr::one())
                        .pi(E::Fr::from(5u64))
                });

                let twenty = composer.arithmetic_gate(|gate| {
                    gate.witness(six, seven, None)
                        .add(E::Fr::one(), E::Fr::one())
                        .fan_in_3(E::Fr::one(), seven)
                });

                // There are quite a few ways to check the equation is correct,
                // depending on your circumstance If we already
                // have the output wire, we can constrain the output of the
                // mul_gate to be equal to it If we do not, we
                // can compute it using an `arithmetic_gate`. If the output
                // is public, we can also constrain the output wire of the mul
                // gate to it. This is what this test does
                let output = composer.arithmetic_gate(|gate| {
                    gate.witness(fourteen, twenty, None).mul(E::Fr::one())
                });

                composer.constrain_to_constant(
                    output,
                    E::Fr::from(280u64),
                    None,
                );
            },
            200,
        );
        assert!(res.is_ok(), "{:?}", res.err().unwrap());
    }

    fn test_correct_add_gate<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let zero = composer.zero_var();
                let one = composer.add_input(E::Fr::one());

                let c = composer.arithmetic_gate(|gate| {
                    gate.witness(one, zero, None)
                        .add(E::Fr::one(), E::Fr::one())
                        .constant(E::Fr::from(2u64))
                });

                composer.constrain_to_constant(c, E::Fr::from(3u64), None);
            },
            32,
        );
        assert!(res.is_ok(), "{:?}", res.err().unwrap());
    }

    fn test_correct_big_add_mul_gate<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                // Verify that (4+5+5) * (6+7+7) + (8*9) = 352
                let four = composer.add_input(E::Fr::from(4u64));
                let five = composer.add_input(E::Fr::from(5u64));
                let six = composer.add_input(E::Fr::from(6u64));
                let seven = composer.add_input(E::Fr::from(7u64));
                let nine = composer.add_input(E::Fr::from(9u64));

                let fourteen = composer.arithmetic_gate(|gate| {
                    gate.witness(four, five, None)
                        .add(E::Fr::one(), E::Fr::one())
                        .fan_in_3(E::Fr::one(), five)
                });

                let twenty = composer.arithmetic_gate(|gate| {
                    gate.witness(six, seven, None)
                        .add(E::Fr::one(), E::Fr::one())
                        .fan_in_3(E::Fr::one(), seven)
                });

                let output = composer.arithmetic_gate(|gate| {
                    gate.witness(fourteen, twenty, None)
                        .mul(E::Fr::one())
                        .fan_in_3(E::Fr::from(8u64), nine)
                });

                composer.constrain_to_constant(
                    output,
                    E::Fr::from(352u64),
                    None,
                );
            },
            200,
        );
        assert!(res.is_ok());
    }

    fn test_correct_big_arith_gate<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                // Verify that (4*5)*6 + 4*7 + 5*8 + 9*10 + 11 = 289
                let a = composer.add_input(E::Fr::from(4u64));
                let b = composer.add_input(E::Fr::from(5u64));
                let q_m = E::Fr::from(6u64);
                let q_l = E::Fr::from(7u64);
                let q_r = E::Fr::from(8u64);
                let d = composer.add_input(E::Fr::from(9u64));
                let q_4 = E::Fr::from(10u64);
                let q_c = E::Fr::from(11u64);

                let output = composer.arithmetic_gate(|gate| {
                    gate.witness(a, b, None)
                        .mul(q_m)
                        .add(q_l, q_r)
                        .fan_in_3(q_4, d)
                        .constant(q_c)
                });

                composer.constrain_to_constant(
                    output,
                    E::Fr::from(289u64),
                    None,
                );
            },
            200,
        );
        assert!(res.is_ok());
    }

    fn test_incorrect_big_arith_gate<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                // Verify that (4*5)*6 + 4*7 + 5*8 + 9*12 + 11 != 289
                let a = composer.add_input(E::Fr::from(4u64));
                let b = composer.add_input(E::Fr::from(5u64));
                let q_m = E::Fr::from(6u64);
                let q_l = E::Fr::from(7u64);
                let q_r = E::Fr::from(8u64);
                let d = composer.add_input(E::Fr::from(9u64));
                let q_4 = E::Fr::from(12u64);
                let q_c = E::Fr::from(11u64);

                let output = composer.arithmetic_gate(|gate| {
                    gate.witness(a, b, None)
                        .mul(q_m)
                        .add(q_l, q_r)
                        .fan_in_3(q_4, d)
                        .constant(q_c)
                });

                composer.constrain_to_constant(
                    output,
                    E::Fr::from(289u64),
                    None,
                );
            },
            200,
        );
        assert!(res.is_err());
    }

    fn test_incorrect_add_mul_gate<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                // Verify that (5+5) * (6+7) != 117
                let five = composer.add_input(E::Fr::from(5u64));
                let six = composer.add_input(E::Fr::from(6u64));
                let seven = composer.add_input(E::Fr::from(7u64));

                let five_plus_five = composer.arithmetic_gate(|gate| {
                    gate.witness(five, five, None)
                        .add(E::Fr::one(), E::Fr::one())
                });

                let six_plus_seven = composer.arithmetic_gate(|gate| {
                    gate.witness(six, seven, None)
                        .add(E::Fr::one(), E::Fr::one())
                });

                let output = composer.arithmetic_gate(|gate| {
                    gate.witness(five_plus_five, six_plus_seven, None)
                        .add(E::Fr::one(), E::Fr::one())
                });

                composer.constrain_to_constant(
                    output,
                    E::Fr::from(117u64),
                    None,
                );
            },
            200,
        );
        assert!(res.is_err());
    }

    // Bls12-381 tests
    batch_test!(
        [
            test_public_inputs,
            test_correct_add_mul_gate,
            test_correct_add_gate,
            test_correct_big_add_mul_gate,
            test_correct_big_arith_gate,
            test_incorrect_add_mul_gate,
            test_incorrect_big_arith_gate
        ],
        [] => (
            Bls12_381,
            ark_ed_on_bls12_381::EdwardsParameters
        )
    );

    // Bls12-377 tests
    batch_test!(
        [
            test_public_inputs,
            test_correct_add_mul_gate,
            test_correct_add_gate,
            test_correct_big_add_mul_gate,
            test_correct_big_arith_gate,
            test_incorrect_add_mul_gate,
            test_incorrect_big_arith_gate
        ],
        [] => (
            Bls12_377,
            ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
