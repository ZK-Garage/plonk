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

use crate::{
    constraint_system::Variable, error::Error, permutation::Permutation,
};

use crate::lookup::LookupTable;
use crate::proof_system::pi::PublicInputs;
use ark_ec::{models::TEModelParameters, ModelParameters};
use ark_ff::{PrimeField, ToConstraintField};
use core::cmp::max;
use core::marker::PhantomData;
use hashbrown::HashMap;
use rand_core::{CryptoRng, RngCore};

/// The StandardComposer is the circuit-builder tool that the `plonk` repository
/// provides to create, stored and transformed circuit descriptions
/// into a [`Proof`](crate::proof_system::Proof) at some point.
///
/// A StandardComposer stores the fullcircuit information, being this one
/// all of the witness and circuit descriptors info (values, positions in the
/// circuits, gates and Wires that occupy..), the public inputs, the connection
/// relationships between the witnesses and how they're represented as Wires (so
/// basically the Permutation argument etc..).
///
/// The StandardComposer also grants us a way to introduce our secret
/// witnesses in the form of a [`Variable`] into the circuit description as well
/// as the public inputs. We can do this with methods like
/// [`StandardComposer::add_input`].
///
/// The StandardComposer also contains as associated functions all the
/// necessary tools to be able to instrument the circuits that the user needs
/// through the addition of gates. There are functions that may add a single
/// arithmetic gate to the circuit [`StandardComposer::arithmetic_gate`] and
/// others that can add several gates to the circuit description such as
/// [`StandardComposer::conditional_select`].
///
/// Each gate or group of gates adds a specific functionality or operation to
/// the circuit description, and so, that's why we can understand
/// the StandardComposer as a builder.
#[derive(derivative::Derivative)]
#[derivative(Debug)]
pub struct StandardComposer<F, P>
where
    F: PrimeField,
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
    // Here we introduce 3 new selectors that will be useful for
    // poseidon hashes.
    /// Selector for for w_l^5
    pub(crate) q_hl: Vec<F>,
    /// Selector for for w_r^5
    pub(crate) q_hr: Vec<F>,
    /// Selector for for w_4^5
    pub(crate) q_h4: Vec<F>,
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
    /// Lookup gate selector
    pub(crate) q_lookup: Vec<F>,

    /// Sparse representation of the Public Inputs linking the positions of the
    /// non-zero ones to it's actual values.
    pub(crate) public_inputs: PublicInputs<F>,
    pub(crate) intended_pi_pos: Vec<usize>,

    // Witness vectors
    /// Left wire witness vector.
    pub(crate) w_l: Vec<Variable>,
    /// Right wire witness vector.
    pub(crate) w_r: Vec<Variable>,
    /// Output wire witness vector.
    pub(crate) w_o: Vec<Variable>,
    /// Fourth wire witness vector.
    pub(crate) w_4: Vec<Variable>,

    /// Public lookup table
    pub(crate) lookup_table: LookupTable<F>,

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
    F: PrimeField,
    P: ModelParameters<BaseField = F>,
{
    /// Returns the length of the circuit that can accommodate the lookup table.
    pub fn total_size(&self) -> usize {
        max(self.n, self.lookup_table.size())
    }

    /// Returns the smallest power of two needed for the circuit.
    pub fn circuit_bound(&self) -> usize {
        self.total_size().next_power_of_two()
    }

    /// Returns a reference to the [`PublicInputs`] stored in the
    /// [`StandardComposer`].
    pub fn get_pi(&self) -> &PublicInputs<F> {
        &self.public_inputs
    }

    /// Insert data in the PI starting at the given position and stores the
    /// occupied positions as intended for public inputs.
    pub(crate) fn add_pi<T>(
        &mut self,
        pos: usize,
        item: &T,
    ) -> Result<(), Error>
    where
        T: ToConstraintField<F>,
    {
        let n_positions = self.public_inputs.add_input(pos, item)?;
        self.intended_pi_pos.extend(pos..(pos + n_positions));
        Ok(())
    }
}

impl<F, P> Default for StandardComposer<F, P>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    #[inline]
    fn default() -> Self {
        Self::new()
    }
}

impl<F, P> StandardComposer<F, P>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
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
            q_hl: Vec::with_capacity(expected_size),
            q_hr: Vec::with_capacity(expected_size),
            q_h4: Vec::with_capacity(expected_size),
            q_arith: Vec::with_capacity(expected_size),
            q_range: Vec::with_capacity(expected_size),
            q_logic: Vec::with_capacity(expected_size),
            q_fixed_group_add: Vec::with_capacity(expected_size),
            q_variable_group_add: Vec::with_capacity(expected_size),
            q_lookup: Vec::with_capacity(expected_size),
            public_inputs: PublicInputs::new(),
            intended_pi_pos: Vec::new(),
            w_l: Vec::with_capacity(expected_size),
            w_r: Vec::with_capacity(expected_size),
            w_o: Vec::with_capacity(expected_size),
            w_4: Vec::with_capacity(expected_size),
            lookup_table: LookupTable::new(),
            zero_var: Variable(0),
            variables: HashMap::with_capacity(expected_size),
            perm: Permutation::new(),
            __: PhantomData::<P>,
        };

        // Reserve the first variable to be zero
        composer.zero_var =
            composer.add_witness_to_circuit_description(F::zero());

        // Add dummy constraints
        composer.add_blinding_factors(&mut rand_core::OsRng);

        composer
    }

    /// Witness representation of zero of the first variable of any circuit
    pub fn zero_var(&self) -> Variable {
        self.zero_var
    }

    /// Add Input first calls the Permutation
    /// to generate and allocate a new [`Variable`] `var`.
    ///
    /// The Composer then links the variable to the [`PrimeField`]
    /// and returns it for its use in the system.
    pub fn add_input(&mut self, s: F) -> Variable {
        // Get a new Variable from the permutation
        let var = self.perm.new_variable();
        // The composer now links the Variable returned from
        // the Permutation to the value F.
        self.variables.insert(var, s);

        var
    }

    /// Adds a width-3 poly gate.
    /// This gate gives total freedom to the end user to implement the
    /// corresponding circuits in the most optimized way possible because
    /// the user has access to the full set of variables, as well as
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
        self.q_lookup.push(F::zero());

        // add high degree selectors
        self.q_hl.push(F::zero());
        self.q_hr.push(F::zero());
        self.q_h4.push(F::zero());

        if let Some(pi) = pi {
            self.add_pi(self.n, &pi).unwrap_or_else(|_| {
                panic!("Could not insert PI {:?} at {}", pi, self.n)
            });
        };

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

    /// A gate which outputs a variable whose value is 1 if
    /// the input is 0 and whose value is 0 otherwise
    pub fn is_zero_with_output(&mut self, a: Variable) -> Variable {
        // Get relevant field values
        let a_value = self.variables.get(&a).unwrap();
        let y_value = a_value.inverse().unwrap_or_else(F::one);

        // This has value 1 if input value is zero, value 0 otherwise
        let b_value = F::one() - *a_value * y_value;

        let y = self.add_input(y_value);
        let b = self.add_input(b_value);
        let zero = self.zero_var();

        // Enforce constraints. The constraint system being used here is
        // a * y + b - 1 = 0
        // a * b = 0
        // where y is auxiliary and b is the boolean (a == 0).
        let _a_times_b = self.arithmetic_gate(|gate| {
            gate.witness(a, b, Some(zero)).mul(F::one())
        });

        let _first_constraint = self.arithmetic_gate(|gate| {
            gate.witness(a, y, Some(zero))
                .mul(F::one())
                .fan_in_3(F::one(), b)
                .constant(-F::one())
        });

        b
    }

    /// A gate which outputs a variable whose value is 1 if the
    /// two input variables have equal values and whose value is 0 otherwise.
    pub fn is_eq_with_output(&mut self, a: Variable, b: Variable) -> Variable {
        let difference = self.arithmetic_gate(|gate| {
            gate.witness(a, b, None).add(F::one(), -F::one())
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
            gate.witness(bit, choice_a, None).mul(F::one())
        });

        // 1 - bit
        let one_min_bit = self.arithmetic_gate(|gate| {
            gate.witness(bit, zero, None)
                .add(-F::one(), F::zero())
                .constant(F::one())
        });

        // (1 - bit) * b
        let one_min_bit_choice_b = self.arithmetic_gate(|gate| {
            gate.witness(one_min_bit, choice_b, None).mul(F::one())
        });

        // [ (1 - bit) * b ] + [ bit * a ]
        self.arithmetic_gate(|gate| {
            gate.witness(one_min_bit_choice_b, bit_times_a, None)
                .add(F::one(), F::one())
        })
    }

    /// Adds the polynomial f(x) = x * a to the circuit description where
    /// `x = bit`. If:
    /// bit == 1 => value,
    /// bit == 0 => 0,
    ///
    /// # Note
    /// The `bit` used as input which is a [`Variable`] should have previously
    /// been constrained to be either 1 or 0 using a bool constrain. See:
    /// [`StandardComposer::boolean_gate`].
    pub fn conditional_select_zero(
        &mut self,
        bit: Variable,
        value: Variable,
    ) -> Variable {
        // returns bit * value
        self.arithmetic_gate(|gate| {
            gate.witness(bit, value, None).mul(F::one())
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

    /// This function adds two dummy gates to the circuit
    /// description which are guaranteed to always satisfy the gate equation.
    /// This function is only used in benchmarking
    pub fn add_dummy_constraints(&mut self) {
        let var_six = self.add_input(F::from(6u64));
        let var_one = self.add_input(F::one());
        let var_seven = self.add_input(F::from(7u64));
        let var_min_twenty = self.add_input(-F::from(20u64));

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
        self.q_lookup.push(F::one());
        // add high degree selectors
        self.q_hl.push(F::zero());
        self.q_hr.push(F::zero());
        self.q_h4.push(F::zero());
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

        self.q_m.push(F::one());
        self.q_l.push(F::one());
        self.q_r.push(F::one());
        self.q_o.push(F::one());
        self.q_c.push(F::from(127u64));
        self.q_4.push(F::zero());
        self.q_arith.push(F::one());
        self.q_range.push(F::zero());
        self.q_logic.push(F::zero());
        self.q_fixed_group_add.push(F::zero());
        self.q_variable_group_add.push(F::zero());
        self.q_lookup.push(F::one());
        // add high degree selectors
        self.q_hl.push(F::zero());
        self.q_hr.push(F::zero());
        self.q_h4.push(F::zero());
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

    /// Adds 3 dummy rows to the lookup table
    /// The first rows match the witness values used for `add_dummy_constraint`
    /// This function is only used for benchmarking
    pub fn add_dummy_lookup_table(&mut self) {
        self.lookup_table.insert_row(
            F::from(6u64),
            F::from(7u64),
            -F::from(20u64),
            F::one(),
        );

        self.lookup_table.insert_row(
            -F::from(20u64),
            F::from(6u64),
            F::from(7u64),
            F::zero(),
        );

        self.lookup_table.insert_row(
            F::from(3u64),
            F::one(),
            F::from(4u64),
            F::from(9u64),
        );
    }

    /// This function is used to add a blinding factors to the witness
    /// and permutation polynomials.
    /// All gate selectors are turned off to guarantee the constraints
    /// are still satisfied.
    pub fn add_blinding_factors<R>(&mut self, rng: &mut R)
    where
        R: CryptoRng + RngCore + ?Sized,
    {
        let mut rand_var_1 = self.zero_var();
        let mut rand_var_2 = self.zero_var();
        // Blinding wires
        for _ in 0..2 {
            rand_var_1 = self.add_input(F::rand(rng));
            rand_var_2 = self.add_input(F::rand(rng));
            let rand_var_3 = self.add_input(F::rand(rng));
            let rand_var_4 = self.add_input(F::rand(rng));

            self.w_l.push(rand_var_1);
            self.w_r.push(rand_var_2);
            self.w_o.push(rand_var_3);
            self.w_4.push(rand_var_4);

            // All selectors fixed to 0 so that the constraints are satisfied
            self.q_m.push(F::zero());
            self.q_l.push(F::zero());
            self.q_r.push(F::zero());
            self.q_o.push(F::zero());
            self.q_c.push(F::zero());
            self.q_4.push(F::zero());
            self.q_arith.push(F::zero());
            self.q_range.push(F::zero());
            self.q_logic.push(F::zero());
            self.q_fixed_group_add.push(F::zero());
            self.q_variable_group_add.push(F::zero());
            self.q_lookup.push(F::zero());
            // add high degree selectors
            self.q_hl.push(F::zero());
            self.q_hr.push(F::zero());
            self.q_h4.push(F::zero());

            self.perm.add_variables_to_map(
                rand_var_1, rand_var_2, rand_var_3, rand_var_4, self.n,
            );
            self.n += 1;
        }

        // Blinding Z
        // We add 2 pairs of equal random points

        self.w_l.push(rand_var_1);
        self.w_r.push(rand_var_2);
        self.w_o.push(self.zero_var());
        self.w_4.push(self.zero_var());

        // All selectors fixed to 0 so that the constraints are satisfied
        self.q_m.push(F::zero());
        self.q_l.push(F::zero());
        self.q_r.push(F::zero());
        self.q_o.push(F::zero());
        self.q_c.push(F::zero());
        self.q_4.push(F::zero());
        self.q_arith.push(F::zero());
        self.q_range.push(F::zero());
        self.q_logic.push(F::zero());
        self.q_fixed_group_add.push(F::zero());
        self.q_variable_group_add.push(F::zero());
        self.q_lookup.push(F::zero());
        // add high degree selectors
        self.q_hl.push(F::zero());
        self.q_hr.push(F::zero());
        self.q_h4.push(F::zero());

        self.perm.add_variables_to_map(
            rand_var_1,
            rand_var_2,
            self.zero_var(),
            self.zero_var(),
            self.n,
        );
        self.n += 1;
    }
    /// Utility function that checks on the "front-end"
    /// side of the PLONK implementation if the identity polynomial
    /// is satisfied for each of the [`StandardComposer`]'s gates.
    ///
    /// The recommended usage is to derive the std output and the std error to a
    /// text file and analyze the gates there.
    ///
    /// # Panic
    /// The function by itself will print each circuit gate info until one of
    /// the gates does not satisfy the equation or there are no more gates. If
    /// the cause is an unsatisfied gate equation, the function will panic.
    #[cfg(feature = "trace")]
    pub fn check_circuit_satisfied(&mut self) {
        use ark_ff::BigInteger;

        use crate::constraint_system::SBOX_ALPHA;
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
        let pi_vec = self.public_inputs.as_evals(self.circuit_bound());
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
            let _qfixed = self.q_fixed_group_add[i];
            let _qvar = self.q_variable_group_add[i];
            let q_hl = self.q_hl[i];
            let q_hr = self.q_hr[i];
            let q_h4 = self.q_h4[i];

            let pi = pi_vec[i];

            let a = w_l[i];
            let a_next = w_l[(i + 1) % self.n];
            let b = w_r[i];
            let b_next = w_r[(i + 1) % self.n];
            let c = w_o[i];
            let d = w_4[i];
            let d_next = w_4[(i + 1) % self.n];

            #[cfg(feature = "trace-print")]
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
            - q_hash_1 -> {:?}\n
            - q_hash_2 -> {:?}\n
            - q_hash_3 -> {:?}\n
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
                q_hl,
                q_hr,
                q_h4,
                qarith,
                qrange,
                qlogic,
                _qfixed,
                _qvar,
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
                    + q_hl * a.pow([SBOX_ALPHA])
                    + q_hr * b.pow([SBOX_ALPHA])
                    + q_h4 * d.pow([SBOX_ALPHA])
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

    /// Get value of a variable that was previously created with
    /// [`add_input`](StandardComposer::add_input),
    /// [`conditional_select`](StandardComposer::conditional_select),
    /// [`is_eq_with_output`](StandardComposer::is_eq_with_output)
    /// or other similar method that returns a `Variable`.
    #[inline]
    pub fn value_of_var(&self, var: Variable) -> F {
        self.variables
            .get(&var)
            .copied()
            .expect("the variable does not exist")
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        batch_test, batch_test_field_params,
        commitment::HomomorphicCommitment,
        constraint_system::helper::*,
        proof_system::{Prover, Verifier},
    };
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use rand_core::OsRng;

    /// Tests that a circuit initially has 3 gates.
    fn test_initial_circuit_size<F, P>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
    {
        // NOTE: Circuit size is n+4 because
        // - We have an extra gate which forces the first witness to be zero.
        //   This is used when the advice wire is not being used.
        // - We have two gates which add random values to blind the wires.
        // - Another gate which adds 2 pairs of equal points to blind the
        //   permutation polynomial
        assert_eq!(4, StandardComposer::<F, P>::new().n)
    }

    /// Tests that an empty circuit proof passes.
    fn test_prove_verify<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        // NOTE: Does nothing except add the dummy constraints.
        let res =
            gadget_tester::<F, P, PC>(|_: &mut StandardComposer<F, P>| {}, 200);
        assert!(res.is_ok());
    }

    fn test_correct_is_zero_with_output<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        // Check that it gives true on zero input:
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                let one = composer.add_input(F::one());
                let is_zero = composer.is_zero_with_output(composer.zero_var());
                composer.assert_equal(is_zero, one);
            },
            32,
        );

        // Check that it gives false on non-zero input:
        let res2 = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                let one = composer.add_input(F::one());
                let is_zero = composer.is_zero_with_output(one);
                composer.assert_equal(is_zero, composer.zero_var());
            },
            32,
        );

        assert!(res.is_ok() && res2.is_ok())
    }

    fn test_correct_is_eq_with_output<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        // Check that it gives true on equal inputs:
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                let one = composer.add_input(F::one());

                let field_element = F::one().double();
                let a = composer.add_input(field_element);
                let b = composer.add_input(field_element);
                let is_eq = composer.is_eq_with_output(a, b);
                composer.assert_equal(is_eq, one);
            },
            32,
        );

        // Check that it gives false on non-equal inputs:
        let res2 = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                let field_element = F::one().double();
                let a = composer.add_input(field_element);
                let b = composer.add_input(field_element.double());
                let is_eq = composer.is_eq_with_output(a, b);
                composer.assert_equal(is_eq, composer.zero_var());
            },
            32,
        );

        assert!(res.is_ok() && res2.is_ok())
    }

    fn test_conditional_select<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
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
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
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

        let public_inputs = prover.cs.get_pi().clone();

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
