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

#[cfg(feature = "trace")]
use ark_ff::BigInteger;
use ark_ff::FftField;
use core::marker::PhantomData;
use hashbrown::HashMap;

use ark_ec::{ModelParameters, SWModelParameters, TEModelParameters};

/// The StandardComposer is the circuit-builder tool that the `dusk-plonk`
/// repository provides so that circuit descriptions can be written, stored and
/// transformed into a [`Proof`](crate::proof_system::Proof) at some point.
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
/// gate to the circuit as for example [`StandardComposer::add_gate`] and others
/// that can add several gates to the circuit description such as
/// [`StandardComposer::conditional_select`].
///
/// Each gate or group of gates adds an specific functionallity or operation to
/// the circuit description, and so, that's why we can understand
/// the StandardComposer as a builder.
#[derive(derivative::Derivative)]
#[derivative(Debug)]
pub struct StandardComposer<F, P>
where
    F: FftField,
    P: ModelParameters<BaseField = F>,
{
    /// Number of arithmetic gates in the circuit
    pub(crate) n: usize,

    // Selector vectors
    /// Multiplier selector
    pub(crate) q_m: Vec<F>,
    /// Left wire selector
    pub(crate) q_l: Vec<F>,
    /// Right wire selector
    pub(crate) q_r: Vec<F>,
    /// Output wire selector
    pub(crate) q_o: Vec<F>,
    /// Fourth wire selector
    pub(crate) q_4: Vec<F>,
    /// Constant wire selector
    pub(crate) q_c: Vec<F>,
    /// Arithmetic wire selector
    pub(crate) q_arith: Vec<F>,
    /// Range selector
    pub(crate) q_range: Vec<F>,
    /// Logic selector
    pub(crate) q_logic: Vec<F>,
    /// Fixed base group addition selector
    pub(crate) q_fixed_group_add: Vec<F>,
    /// Variable base group addition selector
    pub(crate) q_variable_group_add: Vec<F>,

    /// Sparse representation of the Public Inputs linking the positions of the
    /// non-zero ones to it's actual values.
    pub(crate) public_inputs_sparse_store: BTreeMap<usize, F>,

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
    pub(crate) variables: HashMap<Variable, F>,

    /// Permutation argument.
    pub(crate) perm: Permutation,

    /// Type Parameter Marker
    __: PhantomData<P>,
}

impl<F, P> StandardComposer<F, P>
where
    F: FftField,
    P: ModelParameters<BaseField = F>,
{
    /// Returns the number of gates in the circuit
    pub fn circuit_size(&self) -> usize {
        self.n
    }

    /// Constructs a dense vector of the Public Inputs from the positions and
    /// the sparse vector that contains the values.
    pub fn construct_dense_pi_vec(&self) -> Vec<F> {
        let mut pi = vec![F::zero(); self.n];
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

impl<F, P> Default for StandardComposer<F, P>
where
    F: FftField,
    P: ModelParameters<BaseField = F>,
{
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl<F, P> StandardComposer<F, P>
where
    F: FftField,
    P: ModelParameters<BaseField = F>,
{
    /// Generates a new empty `StandardComposer` with all of it's fields
    /// set to hold an initial capacity of 0.
    ///
    /// # Note
    ///
    /// The usage of this may cause lots of re-allocations since the `Composer`
    /// holds `Vec` for every polynomial, and these will need to be re-allocated
    /// each time the circuit grows considerably.
    pub fn new() -> Self {
        Self::with_expected_size(0)
    }

    /// Fixes a [`Variable`] in the witness to be a part of the circuit
    /// description.
    pub fn add_witness_to_circuit_description(&mut self, value: F) -> Variable {
        let var = self.add_input(value);
        self.constrain_to_constant(var, value, None);
        var
    }

    /// Creates a new circuit with an expected circuit size.
    /// This will allow for less reallocations when building the circuit
    /// since the `Vec`s will already have an appropriate allocation at the
    /// beginning of the composing stage.
    pub fn with_expected_size(expected_size: usize) -> Self {
        let mut composer = Self {
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
            __: PhantomData::<P>,
        };

        // Reserve the first variable to be zero
        composer.zero_var =
            composer.add_witness_to_circuit_description(F::zero());

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
    /// The Composer then links the variable to the [`F`]
    /// and returns it for its use in the system.
    ///
    /// [`F`]: Field
    pub fn add_input(&mut self, s: F) -> Variable {
        // Get a new Variable from the permutation
        let var = self.perm.new_variable();
        // The composer now links the F to the Variable returned from
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
        q_m: F,
        q_l: F,
        q_r: F,
        q_o: F,
        q_c: F,
        pi: Option<F>,
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
        self.q_4.push(F::zero());
        self.q_arith.push(F::one());

        self.q_range.push(F::zero());
        self.q_logic.push(F::zero());
        self.q_fixed_group_add.push(F::zero());
        self.q_variable_group_add.push(F::zero());

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
        constant: F,
        pi: Option<F>,
    ) {
        self.poly_gate(
            a,
            a,
            a,
            F::zero(),
            F::one(),
            F::zero(),
            F::zero(),
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
            F::zero(),
            F::one(),
            -F::one(),
            F::zero(),
            F::zero(),
            None,
        );
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
        // bit * choice_a
        let bit_times_a = self.mul(F::one(), bit, choice_a, F::zero(), None);

        // 1 - bit
        let one_min_bit = self.add(
            (-F::one(), bit),
            (F::zero(), self.zero_var),
            F::one(),
            None,
        );

        // (1 - bit) * b
        let one_min_bit_choice_b =
            self.mul(F::one(), one_min_bit, choice_b, F::zero(), None);

        // [ (1 - bit) * b ] + [ bit * a ]
        self.add(
            (F::one(), one_min_bit_choice_b),
            (F::one(), bit_times_a),
            F::zero(),
            None,
        )
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
        self.mul(F::one(), bit, value, F::zero(), None)
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

        let f_x_scalar = F::one() - bit_scalar + (*bit_scalar * value_scalar);
        let f_x = self.add_input(f_x_scalar);

        self.poly_gate(
            bit,
            value,
            f_x,
            F::one(),
            -F::one(),
            F::zero(),
            -F::one(),
            F::one(),
            None,
        );

        f_x
    }

    /// This function is used to add a blinding factor to the witness
    /// polynomials. It essentially adds two dummy gates to the circuit
    /// description which are guaranteed to always satisfy the gate equation.
    pub fn add_dummy_constraints(&mut self) {
        // Add a dummy constraint so that we do not have zero polynomials
        self.q_m.push(F::from(1u64));
        self.q_l.push(F::from(2u64));
        self.q_r.push(F::from(3u64));
        self.q_o.push(F::from(4u64));
        self.q_c.push(F::from(4u64));
        self.q_4.push(F::one());
        self.q_arith.push(F::one());
        self.q_range.push(F::zero());
        self.q_logic.push(F::zero());
        self.q_fixed_group_add.push(F::zero());
        self.q_variable_group_add.push(F::zero());
        let var_six = self.add_input(F::from(6u64));
        let var_one = self.add_input(F::from(1u64));
        let var_seven = self.add_input(F::from(7u64));
        let var_min_twenty = self.add_input(-F::from(20u64));
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
        self.q_m.push(F::from(1u64));
        self.q_l.push(F::from(1u64));
        self.q_r.push(F::from(1u64));
        self.q_o.push(F::from(1u64));
        self.q_c.push(F::from(127u64));
        self.q_4.push(F::zero());
        self.q_arith.push(F::one());
        self.q_range.push(F::zero());
        self.q_logic.push(F::zero());
        self.q_fixed_group_add.push(F::zero());
        self.q_variable_group_add.push(F::zero());
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
        let w_l: Vec<&F> = self
            .w_l
            .iter()
            .map(|w_l_i| self.variables.get(w_l_i).unwrap())
            .collect();
        let w_r: Vec<&F> = self
            .w_r
            .iter()
            .map(|w_r_i| self.variables.get(w_r_i).unwrap())
            .collect();
        let w_o: Vec<&F> = self
            .w_o
            .iter()
            .map(|w_o_i| self.variables.get(w_o_i).unwrap())
            .collect();
        let w_4: Vec<&F> = self
            .w_4
            .iter()
            .map(|w_4_i| self.variables.get(w_4_i).unwrap())
            .collect();
        // Computes f(f-1)(f-2)(f-3)
        let delta = |f: F| -> F {
            let f_1 = f - F::one();
            let f_2 = f - F::from(2u64);
            let f_3 = f - F::from(3u64);
            f * f_1 * f_2 * f_3
        };
        let pi_vec = self.construct_dense_pi_vec();
        let four = F::from(4u64);
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

            let k = qarith
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
                        + match (qlogic == F::one(), qlogic == -F::one()) {
                            (true, false) => {
                                let a_bits = a.into_repr().to_bits_le();
                                let b_bits = b.into_repr().to_bits_le();
                                let a_and_b = a_bits
                                    .iter()
                                    .zip(b_bits)
                                    .map(|(a_bit, b_bit)| a_bit & b_bit)
                                    .collect::<Vec<bool>>();

                                F::from_repr(
                                    <F as PrimeField>::BigInt::from_bits_le(
                                        &a_and_b,
                                    ),
                                )
                                .unwrap()
                                    - *d
                            }
                            (false, true) => {
                                let a_bits = a.into_repr().to_bits_le();
                                let b_bits = b.into_repr().to_bits_le();
                                let a_xor_b = a_bits
                                    .iter()
                                    .zip(b_bits)
                                    .map(|(a_bit, b_bit)| a_bit ^ b_bit)
                                    .collect::<Vec<bool>>();

                                F::from_repr(
                                    <F as PrimeField>::BigInt::from_bits_le(
                                        &a_xor_b,
                                    ),
                                )
                                .unwrap()
                                    - *d
                            }
                            (false, false) => F::zero(),
                            _ => unreachable!(),
                        })
                + qrange
                    * (delta(*c - four * d)
                        + delta(*b - four * c)
                        + delta(*a - four * b)
                        + delta(*d_next - four * a));

            assert_eq!(k, F::zero(), "Check failed at gate {}", i,);
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::commitment::HomomorphicCommitment;
    use crate::constraint_system::helper::*;
    use crate::prelude::Prover;
    use crate::prelude::Verifier;
    use crate::{batch_test, batch_test_field_params};
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_ec::models::TEModelParameters;
    use ark_ec::PairingEngine;
    use ark_ff::{FftField, PrimeField};
    use ark_poly::univariate::DensePolynomial;
    use ark_poly_commit::PolynomialCommitment;
    use rand_core::OsRng;

    /// Tests that a circuit initially has 3 gates.
    fn test_initial_circuit_size<F, P>()
    where
        F: FftField,
        P: ModelParameters<BaseField = F>,
    {
        // NOTE: Circuit size is n+3 because
        // - We have an extra gate which forces the first witness to be zero.
        //   This is used when the advice wire is not being used.
        // - We have two gates which ensure that the permutation polynomial is
        //   not the identity and
        // - Another gate which ensures that the selector polynomials are not
        //   all zeroes
        assert_eq!(3, StandardComposer::<F, P>::new().circuit_size())
    }

    /// Tests that an empty circuit proof passes.
    fn test_prove_verify<F, P, PC>()
    where
        F: FftField + PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: PolynomialCommitment<F, DensePolynomial<F>>
            + HomomorphicCommitment<F>,
    {
        // NOTE: Does nothing except add the dummy constraints.
        let res =
            gadget_tester::<F, P, PC>(|_: &mut StandardComposer<F, P>| {}, 200);
        assert!(res.is_ok());
    }

    fn test_conditional_select<F, P, PC>()
    where
        F: FftField + PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: PolynomialCommitment<F, DensePolynomial<F>>
            + HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                let bit_1 = composer.add_input(F::one());
                let bit_0 = composer.zero_var();

                let choice_a = composer.add_input(F::from(10u64));
                let choice_b = composer.add_input(F::from(20u64));

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
    fn test_multiple_proofs<F, P, PC>()
    where
        F: FftField + PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: PolynomialCommitment<F, DensePolynomial<F>>
            + HomomorphicCommitment<F>,
    {
        let u_params = PC::setup(2 * 30, None, &mut OsRng).unwrap();

        // Create a prover struct
        let mut prover: Prover<F, P, PC> = Prover::new(b"demo");

        // Add gadgets
        dummy_gadget(10, prover.mut_cs());

        // Commit Key
        let (ck, vk) = PC::trim(&u_params, 2 * 20, 0, None).unwrap();

        // Preprocess circuit
        prover.preprocess(&ck).unwrap();

        let public_inputs = prover.cs.construct_dense_pi_vec();

        let mut proofs = Vec::new();

        // Compute multiple proofs
        for _ in 0..3 {
            proofs.push(prover.prove(&ck).unwrap());

            // Add another witness instance
            dummy_gadget(10, prover.mut_cs());
        }

        // Verifier
        //
        let mut verifier = Verifier::<F, P, PC>::new(b"demo");

        // Add gadgets
        dummy_gadget(10, verifier.mut_cs());

        // Preprocess
        verifier.preprocess(&ck).unwrap();

        for proof in proofs {
            assert!(verifier.verify(&proof, &vk, &public_inputs).is_ok());
        }
    }

    // Tests for Bls12_381
    batch_test_field_params!(
        [
            test_initial_circuit_size
        ],
        [] => (
            Bls12_381,
            ark_ed_on_bls12_381::EdwardsParameters

        )
    );

    // Tests for Bls12_377
    batch_test_field_params!(
        [
            test_initial_circuit_size
        ],
        [] => (
            Bls12_377,
            ark_ed_on_bls12_377::EdwardsParameters
        )
    );

    // Tests for Bls12_381
    batch_test!(
        [
            test_prove_verify,
            test_multiple_proofs,
            test_conditional_select ],
        [] => (
            Bls12_381,
            ark_ed_on_bls12_381::EdwardsParameters
        )
    );

    // Tests for Bls12_377
    batch_test!(
        [
            test_prove_verify,
            test_multiple_proofs,
            test_conditional_select
                    ],
        [] => (
            Bls12_377,
            ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
