// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! A `Composer` could be understood as some sort of Trait that is actually
//! defining some kind of Circuit Builder for PLONK.
//!
//! In that sense, here we have the implementation of the [`StandardComposer`]
//! which has been designed in order to provide the maximum amount of
//! performance while having a big scope in utility terms.
//!
//! It allows us not only to build Add and Mul constraints but also to build
//! ECC op. gates, Range checks, Logical gates (Bitwise ops) etc.

use crate::constraint_system::Variable;
use crate::permutation::Permutation;
use alloc::collections::BTreeMap;
use ark_ec::models::TEModelParameters;
use ark_ec::PairingEngine;
use ark_ff::Field;
#[cfg(feature = "trace")]
use ark_ff::{BigInteger, PrimeField};
use core::marker::PhantomData;
use hashbrown::HashMap;
use num_traits::{One, Zero};

/// The StandardComposer is the circuit-builder tool that the `plonk` repository
/// provides so that circuit descriptions can be written, stored and transformed
/// into a [`Proof`](crate::proof_system::Proof) at some point.
///
/// A StandardComposer stores all of the circuit information, being this one
/// all of the witness and circuit descriptors info (values, positions in the
/// circuits, gates and Wires that occupy..), the public inputs, the connection
/// relationships between the witnesses and how they're repesented as Wires (so
/// basically the Permutation argument etc..).
///
/// The StandardComposer also grants us a way to introduce our secret
/// witnesses in a for of a [`Variable`] into the circuit description as well as
/// the public inputs. We can do this with methods like
/// [`StandardComposer::add_input`].
///
/// The StandardComposer also contains as associated functions all the
/// neccessary tools to be able to istrument the circuits that the user needs
/// through the addition of gates. There are functions that may add a single
/// arithmetic gate to the circuit [`StandardComposer::arithmetic_gate`] and
/// others that can add several gates to the circuit description such as
/// [`StandardComposer::conditional_select`].
///
/// Each gate or group of gates adds an specific functionallity or operation to
/// the circuit description, and so, that's why we can understand
/// the StandardComposer as a builder.
#[derive(derivative::Derivative)]
#[derivative(Debug)]
pub struct StandardComposer<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    /// Number of arithmetic gates in the circuit
    pub(crate) n: usize,

    // Selector vectors
    /// Multiplier selector
    pub(crate) q_m: Vec<E::Fr>,
    /// Left wire selector
    pub(crate) q_l: Vec<E::Fr>,
    /// Right wire selector
    pub(crate) q_r: Vec<E::Fr>,
    /// Output wire selector
    pub(crate) q_o: Vec<E::Fr>,
    /// Fourth wire selector
    pub(crate) q_4: Vec<E::Fr>,
    /// Constant wire selector
    pub(crate) q_c: Vec<E::Fr>,
    /// Arithmetic wire selector
    pub(crate) q_arith: Vec<E::Fr>,
    /// Range selector
    pub(crate) q_range: Vec<E::Fr>,
    /// Logic selector
    pub(crate) q_logic: Vec<E::Fr>,
    /// Fixed base group addition selector
    pub(crate) q_fixed_group_add: Vec<E::Fr>,
    /// Variable base group addition selector
    pub(crate) q_variable_group_add: Vec<E::Fr>,

    /// Sparse representation of the Public Inputs linking the positions of the
    /// non-zero ones to it's actual values.
    pub(crate) public_inputs_sparse_store: BTreeMap<usize, E::Fr>,

    // Witness vectors
    /// Left wire witness vector.
    pub(crate) w_l: Vec<Variable>,
    /// Right wire witness vector.
    pub(crate) w_r: Vec<Variable>,
    /// Output wire witness vector.
    pub(crate) w_o: Vec<Variable>,
    /// Fourth wire witness vector.
    pub(crate) w_4: Vec<Variable>,

    /// A zero Variable that is a part of the circuit description.
    /// We reserve a variable to be zero in the system
    /// This is so that when a gate only uses three wires, we set the fourth
    /// wire to be the variable that references zero
    pub(crate) zero_var: Variable,

    /// These are the actual variable values.
    pub(crate) variables: HashMap<Variable, E::Fr>,

    /// Permutation argument.
    pub(crate) perm: Permutation<E::Fr>,

    /// Type Parameter Marker
    __: PhantomData<P>,
}

impl<E, P> StandardComposer<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    /// Returns the number of gates in the circuit
    pub fn circuit_size(&self) -> usize {
        self.n
    }

    /// Constructs a dense vector of the Public Inputs from the positions and
    /// the sparse vector that contains the values.
    pub fn construct_dense_pi_vec(&self) -> Vec<E::Fr> {
        let mut pi = vec![E::Fr::zero(); self.n];
        self.public_inputs_sparse_store
            .iter()
            .for_each(|(pos, value)| {
                pi[*pos] = *value;
            });
        pi
    }

    /// Returns the positions that the Public Inputs occupy in this Composer
    /// instance.
    pub fn pi_positions(&self) -> Vec<usize> {
        // TODO: Find a more performant solution which can return a ref to a Vec
        // or Iterator.
        self.public_inputs_sparse_store.keys().copied().collect()
    }
}

impl<E, P> Default for StandardComposer<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl<E, P> StandardComposer<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    /// Generates a new empty `StandardComposer` with all of it's fields
    /// set to hold an initial capacity of 0.
    ///
    /// # Note
    ///
    /// The usage of this may cause lots of re-allocations since the `Composer`
    /// holds `Vec` for every polynomial, and these will need to be re-allocated
    /// each time the circuit grows considerably.
    pub fn new() -> StandardComposer<E, P> {
        StandardComposer::with_expected_size(0)
    }

    /// Fixes a [`Variable`] in the witness to be a part of the circuit
    /// description.
    pub fn add_witness_to_circuit_description(
        &mut self,
        value: E::Fr,
    ) -> Variable {
        let var = self.add_input(value);
        self.constrain_to_constant(var, value, None);
        var
    }

    /// Creates a new circuit with an expected circuit size.
    /// This will allow for less reallocations when building the circuit
    /// since the `Vec`s will already have an appropriate allocation at the
    /// beginning of the composing stage.
    pub fn with_expected_size(expected_size: usize) -> Self {
        let mut composer = StandardComposer {
            n: 0,
            q_m: Vec::with_capacity(expected_size),
            q_l: Vec::with_capacity(expected_size),
            q_r: Vec::with_capacity(expected_size),
            q_o: Vec::with_capacity(expected_size),
            q_c: Vec::with_capacity(expected_size),
            q_4: Vec::with_capacity(expected_size),
            q_arith: Vec::with_capacity(expected_size),
            q_range: Vec::with_capacity(expected_size),
            q_logic: Vec::with_capacity(expected_size),
            q_fixed_group_add: Vec::with_capacity(expected_size),
            q_variable_group_add: Vec::with_capacity(expected_size),
            public_inputs_sparse_store: BTreeMap::new(),
            w_l: Vec::with_capacity(expected_size),
            w_r: Vec::with_capacity(expected_size),
            w_o: Vec::with_capacity(expected_size),
            w_4: Vec::with_capacity(expected_size),
            zero_var: Variable(0),
            variables: HashMap::with_capacity(expected_size),
            perm: Permutation::new(),
            __: PhantomData,
        };

        // Reserve the first variable to be zero
        composer.zero_var =
            composer.add_witness_to_circuit_description(E::Fr::zero());

        // Add dummy constraints
        composer.add_dummy_constraints();

        composer
    }

    /// Witness representation of zero of the first variable of any circuit
    pub fn zero_var(&self) -> Variable {
        self.zero_var
    }

    /// Add Input first calls the Permutation
    /// to generate and allocate a new [`Variable`] `var`.
    ///
    /// The Composer then links the variable to the [`E::Fr`]
    /// and returns it for its use in the system.
    ///
    /// [`E::Fr`]: PairingEngine::Fr
    pub fn add_input(&mut self, s: E::Fr) -> Variable {
        // Get a new Variable from the permutation
        let var = self.perm.new_variable();
        // The composer now links the E::Fr to the Variable returned from
        // the Permutation
        self.variables.insert(var, s);

        var
    }

    /// Adds a width-3 poly gate.
    /// This gate gives total freedom to the end user to implement the
    /// corresponding circuits in the most optimized way possible because
    /// the under has access to the whole set of variables, as well as
    /// selector coefficients that take part in the computation of the gate
    /// equation.
    ///
    /// The final constraint added will force the following:
    /// `(a * b) * q_m + a * q_l + b * q_r + q_c + PI + q_o * c = 0`.
    pub fn poly_gate(
        &mut self,
        a: Variable,
        b: Variable,
        c: Variable,
        q_m: E::Fr,
        q_l: E::Fr,
        q_r: E::Fr,
        q_o: E::Fr,
        q_c: E::Fr,
        pi: Option<E::Fr>,
    ) -> (Variable, Variable, Variable) {
        self.w_l.push(a);
        self.w_r.push(b);
        self.w_o.push(c);
        self.w_4.push(self.zero_var);
        self.q_l.push(q_l);
        self.q_r.push(q_r);

        // Add selector vectors
        self.q_m.push(q_m);
        self.q_o.push(q_o);
        self.q_c.push(q_c);
        self.q_4.push(E::Fr::zero());
        self.q_arith.push(E::Fr::one());

        self.q_range.push(E::Fr::zero());
        self.q_logic.push(E::Fr::zero());
        self.q_fixed_group_add.push(E::Fr::zero());
        self.q_variable_group_add.push(E::Fr::zero());

        if let Some(pi) = pi {
            assert!(self
                .public_inputs_sparse_store
                .insert(self.n, pi)
                .is_none());
        }

        self.perm
            .add_variables_to_map(a, b, c, self.zero_var, self.n);
        self.n += 1;

        (a, b, c)
    }

    /// Constrain a [`Variable`] to be equal to
    /// a specific constant value which is part of the circuit description and
    /// **NOT** a Public Input. ie. this value will be the same for all of the
    /// circuit instances and [`Proof`](crate::proof_system::Proof)s generated.
    pub fn constrain_to_constant(
        &mut self,
        a: Variable,
        constant: E::Fr,
        pi: Option<E::Fr>,
    ) {
        self.poly_gate(
            a,
            a,
            a,
            E::Fr::zero(),
            E::Fr::one(),
            E::Fr::zero(),
            E::Fr::zero(),
            -constant,
            pi,
        );
    }

    /// Add a constraint into the circuit description that states that two
    /// [`Variable`]s are equal.
    pub fn assert_equal(&mut self, a: Variable, b: Variable) {
        self.poly_gate(
            a,
            b,
            self.zero_var,
            E::Fr::zero(),
            E::Fr::one(),
            -E::Fr::one(),
            E::Fr::zero(),
            E::Fr::zero(),
            None,
        );
    }

    /// A gate which outputs a variable whose value is 1 if
    /// the input is 0 and whose value is 0 otherwise
    pub fn is_zero_with_output(&mut self, a: Variable) -> Variable {
        // Get relevant field values
        let a_value = self.variables.get(&a).unwrap();
        let y_value = a_value.inverse().unwrap_or_else(E::Fr::one);

        // This has value 1 if input value is zero, value 0 otherwise
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

        let _first_constraint = self.arithmetic_gate(|gate| {
            gate.witness(a, y, Some(zero))
                .mul(E::Fr::one())
                .fan_in_3(E::Fr::one(), b)
                .constant(-E::Fr::one())
        });

        b
    }

    /// A gate which outputs a variable whose value is 1 if the
    /// two input variables have equal values and whose value is 0 otherwise.
    pub fn is_eq_with_output(&mut self, a: Variable, b: Variable) -> Variable {
        let difference = self.arithmetic_gate(|gate| {
            gate.witness(a, b, None).add(E::Fr::one(), -E::Fr::one())
        });
        self.is_zero_with_output(difference)
    }

    /// Conditionally selects a [`Variable`] based on an input bit.
    ///
    /// If:
    /// bit == 1 => choice_a,
    /// bit == 0 => choice_b,
    ///
    /// # Note
    /// The `bit` used as input which is a [`Variable`] should had previously
    /// been constrained to be either 1 or 0 using a bool constrain. See:
    /// [`StandardComposer::boolean_gate`].
    pub fn conditional_select(
        &mut self,
        bit: Variable,
        choice_a: Variable,
        choice_b: Variable,
    ) -> Variable {
        let zero = self.zero_var;
        // bit * choice_a
        let bit_times_a = self.arithmetic_gate(|gate| {
            gate.witness(bit, choice_a, None).mul(E::Fr::one())
        });

        // 1 - bit
        let one_min_bit = self.arithmetic_gate(|gate| {
            gate.witness(bit, zero, None)
                .add(-E::Fr::one(), E::Fr::zero())
                .constant(E::Fr::one())
        });

        // (1 - bit) * b
        let one_min_bit_choice_b = self.arithmetic_gate(|gate| {
            gate.witness(one_min_bit, choice_b, None).mul(E::Fr::one())
        });

        // [ (1 - bit) * b ] + [ bit * a ]
        self.arithmetic_gate(|gate| {
            gate.witness(one_min_bit_choice_b, bit_times_a, None)
                .add(E::Fr::one(), E::Fr::one())
        })
    }

    /// Adds the polynomial f(x) = x * a to the circuit description where
    /// `x = bit`. If:
    /// bit == 1 => value,
    /// bit == 0 => 0,
    ///
    /// # Note
    /// The `bit` used as input which is a [`Variable`] should had previously
    /// been constrained to be either 1 or 0 using a bool constrain. See:
    /// [`StandardComposer::boolean_gate`].
    pub fn conditional_select_zero(
        &mut self,
        bit: Variable,
        value: Variable,
    ) -> Variable {
        // returns bit * value
        self.arithmetic_gate(|gate| {
            gate.witness(bit, value, None).mul(E::Fr::one())
        })
    }

    /// Adds the polynomial f(x) = 1 - x + xa to the circuit description where
    /// `x = bit`. If:
    /// bit == 1 => value,
    /// bit == 0 => 1,
    ///
    /// # Note
    /// The `bit` used as input which is a [`Variable`] should had previously
    /// been constrained to be either 1 or 0 using a bool constrain. See:
    /// [`StandardComposer::boolean_gate`].
    pub fn conditional_select_one(
        &mut self,
        bit: Variable,
        value: Variable,
    ) -> Variable {
        let value_scalar = self.variables.get(&value).unwrap();
        let bit_scalar = self.variables.get(&bit).unwrap();

        let f_x_scalar =
            E::Fr::one() - bit_scalar + (*bit_scalar * value_scalar);
        let f_x = self.add_input(f_x_scalar);

        self.poly_gate(
            bit,
            value,
            f_x,
            E::Fr::one(),
            -E::Fr::one(),
            E::Fr::zero(),
            -E::Fr::one(),
            E::Fr::one(),
            None,
        );

        f_x
    }

    /// This function is used to add a blinding factor to the witness
    /// polynomials. It essentially adds two dummy gates to the circuit
    /// description which are guaranteed to always satisfy the gate equation.
    pub fn add_dummy_constraints(&mut self) {
        // Add a dummy constraint so that we do not have zero polynomials
        self.q_m.push(E::Fr::from(1u64));
        self.q_l.push(E::Fr::from(2u64));
        self.q_r.push(E::Fr::from(3u64));
        self.q_o.push(E::Fr::from(4u64));
        self.q_c.push(E::Fr::from(4u64));
        self.q_4.push(E::Fr::one());
        self.q_arith.push(E::Fr::one());
        self.q_range.push(E::Fr::zero());
        self.q_logic.push(E::Fr::zero());
        self.q_fixed_group_add.push(E::Fr::zero());
        self.q_variable_group_add.push(E::Fr::zero());
        let var_six = self.add_input(E::Fr::from(6u64));
        let var_one = self.add_input(E::Fr::from(1u64));
        let var_seven = self.add_input(E::Fr::from(7u64));
        let var_min_twenty = self.add_input(-E::Fr::from(20u64));
        self.w_l.push(var_six);
        self.w_r.push(var_seven);
        self.w_o.push(var_min_twenty);
        self.w_4.push(var_one);
        self.perm.add_variables_to_map(
            var_six,
            var_seven,
            var_min_twenty,
            var_one,
            self.n,
        );
        self.n += 1;
        //Add another dummy constraint so that we do not get the identity
        // permutation
        self.q_m.push(E::Fr::from(1u64));
        self.q_l.push(E::Fr::from(1u64));
        self.q_r.push(E::Fr::from(1u64));
        self.q_o.push(E::Fr::from(1u64));
        self.q_c.push(E::Fr::from(127u64));
        self.q_4.push(E::Fr::zero());
        self.q_arith.push(E::Fr::one());
        self.q_range.push(E::Fr::zero());
        self.q_logic.push(E::Fr::zero());
        self.q_fixed_group_add.push(E::Fr::zero());
        self.q_variable_group_add.push(E::Fr::zero());
        self.w_l.push(var_min_twenty);
        self.w_r.push(var_six);
        self.w_o.push(var_seven);
        self.w_4.push(self.zero_var);
        self.perm.add_variables_to_map(
            var_min_twenty,
            var_six,
            var_seven,
            self.zero_var,
            self.n,
        );
        self.n += 1;
    }

    /// Utility function that allows to check on the "front-end"
    /// side of the PLONK implementation if the identity polynomial
    /// is satisfied for each one of the [`StandardComposer`]'s gates.
    ///
    /// The recommended usage is to derive the std output and the std error to a
    /// text file and analyze there the gates.
    ///
    /// # Panic
    /// The function by itself will print each circuit gate info until one of
    /// the gates does not satisfy the equation or there are no more gates. If
    /// the cause is an unsatisfied gate equation, the function will panic.
    #[cfg(feature = "trace")]
    pub fn check_circuit_satisfied(&self) {
        let w_l: Vec<&E::Fr> = self
            .w_l
            .iter()
            .map(|w_l_i| self.variables.get(w_l_i).unwrap())
            .collect();
        let w_r: Vec<&E::Fr> = self
            .w_r
            .iter()
            .map(|w_r_i| self.variables.get(w_r_i).unwrap())
            .collect();
        let w_o: Vec<&E::Fr> = self
            .w_o
            .iter()
            .map(|w_o_i| self.variables.get(w_o_i).unwrap())
            .collect();
        let w_4: Vec<&E::Fr> = self
            .w_4
            .iter()
            .map(|w_4_i| self.variables.get(w_4_i).unwrap())
            .collect();
        // Computes f(f-1)(f-2)(f-3)
        let delta = |f: E::Fr| -> E::Fr {
            let f_1 = f - E::Fr::one();
            let f_2 = f - E::Fr::from(2u64);
            let f_3 = f - E::Fr::from(3u64);
            f * f_1 * f_2 * f_3
        };
        let pi_vec = self.construct_dense_pi_vec();
        let four = E::Fr::from(4u64);
        for i in 0..self.n {
            let qm = self.q_m[i];
            let ql = self.q_l[i];
            let qr = self.q_r[i];
            let qo = self.q_o[i];
            let qc = self.q_c[i];
            let q4 = self.q_4[i];
            let qarith = self.q_arith[i];
            let qrange = self.q_range[i];
            let qlogic = self.q_logic[i];
            #[cfg(all(feature = "trace-print", feature = "std"))]
            let qfixed = self.q_fixed_group_add[i];
            #[cfg(all(feature = "trace-print", feature = "std"))]
            let qvar = self.q_variable_group_add[i];
            let pi = pi_vec[i];

            let a = w_l[i];
            let a_next = w_l[(i + 1) % self.n];
            let b = w_r[i];
            let b_next = w_r[(i + 1) % self.n];
            let c = w_o[i];
            let d = w_4[i];
            let d_next = w_4[(i + 1) % self.n];

            #[cfg(all(feature = "trace-print", feature = "std"))]
            std::println!(
                "--------------------------------------------\n
            #Gate Index = {}
            #Selector Polynomials:\n
            - qm -> {:?}\n
            - ql -> {:?}\n
            - qr -> {:?}\n
            - q4 -> {:?}\n
            - qo -> {:?}\n
            - qc -> {:?}\n
            - q_arith -> {:?}\n
            - q_range -> {:?}\n
            - q_logic -> {:?}\n
            - q_fixed_group_add -> {:?}\n
            - q_variable_group_add -> {:?}\n
            # Witness polynomials:\n
            - w_l -> {:?}\n
            - w_r -> {:?}\n
            - w_o -> {:?}\n
            - w_4 -> {:?}\n",
                i,
                qm,
                ql,
                qr,
                q4,
                qo,
                qc,
                qarith,
                qrange,
                qlogic,
                qfixed,
                qvar,
                a,
                b,
                c,
                d
            );

            let k =
                qarith
                    * ((qm * a * b)
                        + (ql * a)
                        + (qr * b)
                        + (qo * c)
                        + (q4 * d)
                        + pi
                        + qc)
                    + qlogic
                        * (((delta(*a_next - four * a)
                            - delta(*b_next - four * b))
                            * c)
                            + delta(*a_next - four * a)
                            + delta(*b_next - four * b)
                            + delta(*d_next - four * d)
                            + match (
                                qlogic == E::Fr::one(),
                                qlogic == -E::Fr::one(),
                            ) {
                                (true, false) => {
                                    let a_bits = a.into_repr().to_bits_le();
                                    let b_bits = b.into_repr().to_bits_le();
                                    let a_and_b = a_bits
                                        .iter()
                                        .zip(b_bits)
                                        .map(|(a_bit, b_bit)| a_bit & b_bit)
                                        .collect::<Vec<bool>>();

                                    E::Fr::from_repr(
                                    <E::Fr as PrimeField>::BigInt::from_bits_le(
                                        &a_and_b,
                                    ),
                                ).unwrap() - *d
                                }
                                (false, true) => {
                                    let a_bits = a.into_repr().to_bits_le();
                                    let b_bits = b.into_repr().to_bits_le();
                                    let a_xor_b = a_bits
                                        .iter()
                                        .zip(b_bits)
                                        .map(|(a_bit, b_bit)| a_bit ^ b_bit)
                                        .collect::<Vec<bool>>();

                                    E::Fr::from_repr(
                                    <E::Fr as PrimeField>::BigInt::from_bits_le(
                                        &a_xor_b,
                                    ),
                                ).unwrap() - *d
                                }
                                (false, false) => E::Fr::zero(),
                                _ => unreachable!(),
                            })
                    + qrange
                        * (delta(*c - four * d)
                            + delta(*b - four * c)
                            + delta(*a - four * b)
                            + delta(*d_next - four * a));

            assert_eq!(k, E::Fr::zero(), "Check failed at gate {}", i,);
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::batch_test;
    use crate::constraint_system::helper::*;
    use crate::proof_system::{Prover, Verifier};
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly_commit::kzg10::{self, Powers, UniversalParams, KZG10};
    use ark_poly_commit::sonic_pc::SonicKZG10;
    use ark_poly_commit::PolynomialCommitment;
    use rand_core::OsRng;

    /// Tests that a circuit initially has 3 gates.
    fn test_initial_circuit_size<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        // NOTE: Circuit size is n+3 because
        // - We have an extra gate which forces the first witness to be zero.
        //   This is used when the advice wire is not being used.
        // - We have two gates which ensure that the permutation polynomial is
        //   not the identity and
        // - Another gate which ensures that the selector polynomials are not
        //   all zeroes
        assert_eq!(3, StandardComposer::<E, P>::new().circuit_size())
    }

    /// Tests that an empty circuit proof passes.
    fn test_prove_verify<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        // NOTE: Does nothing except add the dummy constraints.
        let res = gadget_tester(|_: &mut StandardComposer<E, P>| {}, 200);
        assert!(res.is_ok());
    }

    fn test_correct_is_zero_with_output<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        // Check that it gives true on zero input:
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let one = composer.add_input(E::Fr::one());
                let is_zero = composer.is_zero_with_output(composer.zero_var());
                composer.assert_equal(is_zero, one);
            },
            32,
        );

        // Check that it gives false on non-zero input:
        let res2 = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let one = composer.add_input(E::Fr::one());
                let is_zero = composer.is_zero_with_output(one);
                composer.assert_equal(is_zero, composer.zero_var());
            },
            32,
        );

        assert!(res.is_ok() && res2.is_ok())
    }

    fn test_correct_is_eq_with_output<E, P>()
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
                let is_eq = composer.is_eq_with_output(a, b);
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
                let is_eq = composer.is_eq_with_output(a, b);
                composer.assert_equal(is_eq, composer.zero_var());
            },
            32,
        );

        assert!(res.is_ok() && res2.is_ok())
    }

    fn test_conditional_select<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let bit_1 = composer.add_input(E::Fr::one());
                let bit_0 = composer.zero_var();

                let choice_a = composer.add_input(E::Fr::from(10u64));
                let choice_b = composer.add_input(E::Fr::from(20u64));

                let choice =
                    composer.conditional_select(bit_1, choice_a, choice_b);
                composer.assert_equal(choice, choice_a);

                let choice =
                    composer.conditional_select(bit_0, choice_a, choice_b);
                composer.assert_equal(choice, choice_b);
            },
            32,
        );
        assert!(res.is_ok(), "{:?}", res.err().unwrap());
    }

    // FIXME: Move this to integration tests
    fn test_multiple_proofs<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let u_params: UniversalParams<E> =
            KZG10::<E, DensePolynomial<E::Fr>>::setup(
                2 * 30,
                false,
                &mut OsRng,
            )
            .unwrap();

        // Create a prover struct
        let mut prover: Prover<E, P> = Prover::new(b"demo");

        // Add gadgets
        dummy_gadget(10, prover.mut_cs());

        // Commit Key
        let (ck, _) = SonicKZG10::<E, DensePolynomial<E::Fr>>::trim(
            &u_params,
            2 * 20,
            0,
            None,
        )
        .unwrap();
        let powers = Powers {
            powers_of_g: ck.powers_of_g.into(),
            powers_of_gamma_g: ck.powers_of_gamma_g.into(),
        };

        // Preprocess circuit
        prover.preprocess(&powers).unwrap();

        let public_inputs = prover.cs.construct_dense_pi_vec();

        let mut proofs = Vec::new();

        // Compute multiple proofs
        for _ in 0..3 {
            proofs.push(prover.prove(&powers).unwrap());

            // Add another witness instance
            dummy_gadget(10, prover.mut_cs());
        }

        // Verifier
        //
        let mut verifier: Verifier<E, P> = Verifier::new(b"demo");

        // Add gadgets
        dummy_gadget(10, verifier.mut_cs());

        // Commit and Verifier Key
        let (sonic_ck, sonic_vk) =
            SonicKZG10::<E, DensePolynomial<E::Fr>>::trim(
                &u_params,
                2 * 20,
                0,
                None,
            )
            .unwrap();
        let powers = Powers {
            powers_of_g: sonic_ck.powers_of_g.into(),
            powers_of_gamma_g: sonic_ck.powers_of_gamma_g.into(),
        };

        let vk = kzg10::VerifierKey {
            g: sonic_vk.g,
            gamma_g: sonic_vk.gamma_g,
            h: sonic_vk.h,
            beta_h: sonic_vk.beta_h,
            prepared_h: sonic_vk.prepared_h,
            prepared_beta_h: sonic_vk.prepared_beta_h,
        };

        // Preprocess
        verifier.preprocess(&powers).unwrap();

        for proof in proofs {
            assert!(verifier.verify(&proof, &vk, &public_inputs).is_ok());
        }
    }

    // Tests for Bls12_381
    batch_test!(
        [
            test_initial_circuit_size,
            test_prove_verify,
            test_correct_is_zero_with_output,
            test_correct_is_eq_with_output,
            test_conditional_select,
            test_multiple_proofs
        ],
        [] => (
            Bls12_381,
            ark_ed_on_bls12_381::EdwardsParameters
        )
    );

    // Tests for Bls12_377
    batch_test!(
        [
            test_initial_circuit_size,
            test_prove_verify,
            test_correct_is_zero_with_output,
            test_correct_is_eq_with_output,
            test_conditional_select,
            test_multiple_proofs
        ],
        [] => (
            Bls12_377,
            ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
