// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Simple Arithmetic Gates

use crate::constraint_system::{StandardComposer, Variable};
use ark_ec::TEModelParameters;
use ark_ff::PrimeField;

#[derive(Debug, Clone, Copy)]
pub struct ArithmeticGate<F>
where
    F: PrimeField,
{
    pub(crate) witness: Option<(Variable, Variable, Option<Variable>)>,
    pub(crate) fan_in_3: Option<(F, Variable)>,
    pub(crate) mul_selector: F,
    pub(crate) add_selectors: (F, F),
    pub(crate) out_selector: F,
    pub(crate) const_selector: F,
    pub(crate) pi: Option<F>,
}

impl<F> Default for ArithmeticGate<F>
where
    F: PrimeField,
{
    fn default() -> Self {
        Self {
            witness: None,
            fan_in_3: None,
            mul_selector: F::zero(),
            add_selectors: (F::zero(), F::zero()),
            out_selector: -F::one(),
            const_selector: F::zero(),
            pi: None,
        }
    }
}

impl<F> ArithmeticGate<F>
where
    F: PrimeField,
{
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

    pub fn fan_in_3(&mut self, q_4: F, w_4: Variable) -> &mut Self {
        self.fan_in_3 = Some((q_4, w_4));
        self
    }

    pub fn mul(&mut self, q_m: F) -> &mut Self {
        self.mul_selector = q_m;
        self
    }

    pub fn add(&mut self, q_l: F, q_r: F) -> &mut Self {
        self.add_selectors = (q_l, q_r);
        self
    }

    pub fn out(&mut self, q_o: F) -> &mut Self {
        self.out_selector = q_o;
        self
    }

    pub fn constant(&mut self, q_c: F) -> &mut Self {
        self.const_selector = q_c;
        self
    }

    pub fn pi(&mut self, pi: F) -> &mut Self {
        self.pi = Some(pi);
        self
    }

    pub fn build(&mut self) -> Self {
        *self
    }
}

impl<F, P> StandardComposer<F, P>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    /// Function used to generate any arithmetic gate with fan-in-2 or fan-in-3.
    pub fn arithmetic_gate<Fn>(&mut self, func: Fn) -> Variable
    where
        Fn: FnOnce(&mut ArithmeticGate<F>) -> &mut ArithmeticGate<F>,
    {
        let gate = {
            let mut gate = ArithmeticGate::<F>::new();
            func(&mut gate).build()
        };

        if gate.witness.is_none() {
            panic!("Missing left and right wire witnesses")
        }

        let (q4, w4) = gate.fan_in_3.unwrap_or((F::zero(), self.zero_var));
        self.w_4.push(w4);
        self.q_4.push(q4);

        let gate_witness = gate.witness.unwrap();
        self.w_l.push(gate_witness.0);
        self.w_r.push(gate_witness.1);
        self.q_l.push(gate.add_selectors.0);
        self.q_r.push(gate.add_selectors.1);

        // Add selector vectors
        self.q_m.push(gate.mul_selector);
        self.q_o.push(gate.out_selector);
        self.q_c.push(gate.const_selector);

        // add high degree selectors
        self.q_hl.push(F::zero());
        self.q_hr.push(F::zero());
        self.q_h4.push(F::zero());

        self.q_arith.push(F::one());
        self.q_range.push(F::zero());
        self.q_logic.push(F::zero());
        self.q_fixed_group_add.push(F::zero());
        self.q_variable_group_add.push(F::zero());
        self.q_lookup.push(F::zero());

        if let Some(pi) = gate.pi {
            self.add_pi(self.n, &pi).unwrap_or_else(|_| {
                panic!("Could not insert PI {:?} at {}", pi, self.n)
            });
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
    use crate::{
        batch_test, commitment::HomomorphicCommitment,
        constraint_system::helper::*,
    };
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;

    fn test_public_inputs<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                let var_one = composer.add_input(F::one());

                let should_be_three = composer.arithmetic_gate(|gate| {
                    gate.witness(var_one, var_one, None)
                        .add(F::one(), F::one())
                        .pi(F::one())
                });

                composer.constrain_to_constant(
                    should_be_three,
                    F::from(3u64),
                    None,
                );

                let should_be_four = composer.arithmetic_gate(|gate| {
                    gate.witness(var_one, var_one, None)
                        .add(F::one(), F::one())
                        .pi(F::from(2u64))
                });

                composer.constrain_to_constant(
                    should_be_four,
                    F::from(4u64),
                    None,
                );
            },
            200,
        );
        assert!(res.is_ok(), "{:?}", res.err().unwrap());
    }

    fn test_correct_add_mul_gate<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                // Verify that (4+5+5) * (6+7+7) = 280
                let four = composer.add_input(F::from(4u64));
                let five = composer.add_input(F::from(5u64));
                let six = composer.add_input(F::from(6u64));
                let seven = composer.add_input(F::from(7u64));

                let fourteen = composer.arithmetic_gate(|gate| {
                    gate.witness(four, five, None)
                        .add(F::one(), F::one())
                        .pi(F::from(5u64))
                });

                let twenty = composer.arithmetic_gate(|gate| {
                    gate.witness(six, seven, None)
                        .add(F::one(), F::one())
                        .fan_in_3(F::one(), seven)
                });

                // There are quite a few ways to check the equation is correct,
                // depending on your circumstance If we already
                // have the output wire, we can constrain the output of the
                // mul_gate to be equal to it If we do not, we
                // can compute it using an `arithmetic_gate`. If the output
                // is public, we can also constrain the output wire of the mul
                // gate to it. This is what this test does
                let output = composer.arithmetic_gate(|gate| {
                    gate.witness(fourteen, twenty, None).mul(F::one())
                });

                composer.constrain_to_constant(output, F::from(280u64), None);
            },
            200,
        );
        assert!(res.is_ok(), "{:?}", res.err().unwrap());
    }

    fn test_correct_add_gate<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                let zero = composer.zero_var();
                let one = composer.add_input(F::one());

                let c = composer.arithmetic_gate(|gate| {
                    gate.witness(one, zero, None)
                        .add(F::one(), F::one())
                        .constant(F::from(2u64))
                });

                composer.constrain_to_constant(c, F::from(3u64), None);
            },
            32,
        );
        assert!(res.is_ok(), "{:?}", res.err().unwrap());
    }

    fn test_correct_big_add_mul_gate<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                // Verify that (4+5+5) * (6+7+7) + (8*9) = 352
                let four = composer.add_input(F::from(4u64));
                let five = composer.add_input(F::from(5u64));
                let six = composer.add_input(F::from(6u64));
                let seven = composer.add_input(F::from(7u64));
                let nine = composer.add_input(F::from(9u64));

                let fourteen = composer.arithmetic_gate(|gate| {
                    gate.witness(four, five, None)
                        .add(F::one(), F::one())
                        .fan_in_3(F::one(), five)
                });

                let twenty = composer.arithmetic_gate(|gate| {
                    gate.witness(six, seven, None)
                        .add(F::one(), F::one())
                        .fan_in_3(F::one(), seven)
                });

                let output = composer.arithmetic_gate(|gate| {
                    gate.witness(fourteen, twenty, None)
                        .mul(F::one())
                        .fan_in_3(F::from(8u64), nine)
                });

                composer.constrain_to_constant(output, F::from(352u64), None);
            },
            200,
        );
        assert!(res.is_ok());
    }

    fn test_correct_big_arith_gate<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                // Verify that (4*5)*6 + 4*7 + 5*8 + 9*10 + 11 = 289
                let a = composer.add_input(F::from(4u64));
                let b = composer.add_input(F::from(5u64));
                let q_m = F::from(6u64);
                let q_l = F::from(7u64);
                let q_r = F::from(8u64);
                let d = composer.add_input(F::from(9u64));
                let q_4 = F::from(10u64);
                let q_c = F::from(11u64);

                let output = composer.arithmetic_gate(|gate| {
                    gate.witness(a, b, None)
                        .mul(q_m)
                        .add(q_l, q_r)
                        .fan_in_3(q_4, d)
                        .constant(q_c)
                });

                composer.constrain_to_constant(output, F::from(289u64), None);
            },
            200,
        );
        assert!(res.is_ok());
    }

    fn test_incorrect_big_arith_gate<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                // Verify that (4*5)*6 + 4*7 + 5*8 + 9*12 + 11 != 289
                let a = composer.add_input(F::from(4u64));
                let b = composer.add_input(F::from(5u64));
                let q_m = F::from(6u64);
                let q_l = F::from(7u64);
                let q_r = F::from(8u64);
                let d = composer.add_input(F::from(9u64));
                let q_4 = F::from(12u64);
                let q_c = F::from(11u64);

                let output = composer.arithmetic_gate(|gate| {
                    gate.witness(a, b, None)
                        .mul(q_m)
                        .add(q_l, q_r)
                        .fan_in_3(q_4, d)
                        .constant(q_c)
                });

                composer.constrain_to_constant(output, F::from(289u64), None);
            },
            200,
        );
        assert!(res.is_err());
    }

    fn test_incorrect_add_mul_gate<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                // Verify that (5+5) * (6+7) != 117
                let five = composer.add_input(F::from(5u64));
                let six = composer.add_input(F::from(6u64));
                let seven = composer.add_input(F::from(7u64));

                let five_plus_five = composer.arithmetic_gate(|gate| {
                    gate.witness(five, five, None).add(F::one(), F::one())
                });

                let six_plus_seven = composer.arithmetic_gate(|gate| {
                    gate.witness(six, seven, None).add(F::one(), F::one())
                });

                let output = composer.arithmetic_gate(|gate| {
                    gate.witness(five_plus_five, six_plus_seven, None)
                        .add(F::one(), F::one())
                });

                composer.constrain_to_constant(output, F::from(117u64), None);
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
            Bls12_381, ark_ed_on_bls12_381::EdwardsParameters
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
            Bls12_377, ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
