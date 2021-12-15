// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-INFRA. All rights reserved.

//! Variable-base Scalar Multiplication Gate

use crate::constraint_system::ecc::Point;
use crate::constraint_system::{variable::Variable, StandardComposer};
use ark_ec::models::TEModelParameters;
use ark_ec::PairingEngine;
use ark_ff::{BigInteger, Field, FpParameters, PrimeField};
use num_traits::{One, Zero};

impl<E, P> StandardComposer<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    /// Adds a variable-base scalar multiplication to the circuit description.
    ///
    /// # Note
    ///
    /// If you're planning to multiply always by the generator of the scalar
    /// field, you should use [`StandardComposer::fixed_base_scalar_mul`]
    /// which is optimized for fixed_base ops.
    pub fn variable_base_scalar_mul(
        &mut self,
        curve_var: Variable,
        point: Point<E, P>,
    ) -> Point<E, P> {
        // Turn scalar into bits
        let raw_scalar = *self
            .variables
            .get(&curve_var)
            // We can unwrap safely here since it should be impossible to obtain
            // a `Variable` without first linking it inside of the
            // HashMap from which we are calling the `get()` now. Therefore, if
            // the `get()` fn fails now, somethig is going really
            // bad.
            .expect("Variable in existance without referenced scalar");
        let scalar_bits_var = self.scalar_decomposition(curve_var, raw_scalar);

        let identity = Point::identity(self);
        let mut result = identity;

        for bit in scalar_bits_var.into_iter().rev() {
            result = self.point_addition_gate(result, result);
            let point_to_add = self.conditional_select_identity(bit, point);
            result = self.point_addition_gate(result, point_to_add);
        }

        result
    }

    fn scalar_decomposition(
        &mut self,
        witness_var: Variable,
        witness_scalar: E::Fr,
    ) -> Vec<Variable> {
        // Decompose the bits
        let scalar_bits_iter = witness_scalar.into_repr().to_bits_le();

        // Add all the bits into the composer
        let scalar_bits_var: Vec<Variable> = scalar_bits_iter
            .iter()
            .map(|bit| self.add_input(E::Fr::from(*bit as u64)))
            .collect();

        // Take the first 252 bits
        let scalar_bits_var = scalar_bits_var
            [..<P::BaseField as PrimeField>::Params::MODULUS_BITS as usize]
            .to_vec();

        // Now ensure that the bits correctly accumulate to the witness given
        let mut accumulator_var = self.zero_var;
        let mut accumulator_scalar = E::Fr::zero();

        for (power, bit) in scalar_bits_var.iter().enumerate() {
            self.boolean_gate(*bit);

            let two_pow = E::Fr::from(2u64).pow([power as u64, 0, 0, 0]);

            accumulator_var = self.arithmetic_gate(|gate| {
                gate.witness((*bit, accumulator_var, None))
                    .add((two_pow, E::Fr::one()))
            });

            accumulator_scalar +=
                two_pow * E::Fr::from(scalar_bits_iter[power] as u64);
        }
        self.assert_equal(accumulator_var, witness_var);

        scalar_bits_var
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{batch_test, constraint_system::helper::*, util};
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_ec::{twisted_edwards_extended::GroupAffine, AffineCurve};

    fn test_var_base_scalar_mul<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let scalar = E::Fr::from_le_bytes_mod_order(&[
                    182, 44, 247, 214, 94, 14, 151, 208, 130, 16, 200, 204,
                    147, 32, 104, 166, 0, 59, 52, 1, 1, 59, 103, 6, 169, 175,
                    51, 101, 234, 180, 125, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0,
                ]);
                let secret_scalar = composer.add_input(scalar);

                let (x, y) = P::AFFINE_GENERATOR_COEFFS;
                let generator = GroupAffine::new(x, y);

                let expected_point: GroupAffine<P> = AffineCurve::mul(
                    &generator,
                    util::to_embedded_curve_scalar::<E, P>(scalar),
                )
                .into();

                let point = composer.add_affine(generator);

                let point_scalar =
                    composer.variable_base_scalar_mul(secret_scalar, point);

                composer
                    .assert_equal_public_point(point_scalar, expected_point);
            },
            4096,
        );
        assert!(res.is_ok());
    }

    // Tests for Bls12_381
    batch_test!(
        [test_var_base_scalar_mul],
        [] => (
            Bls12_381,
            ark_ed_on_bls12_381::EdwardsParameters
        )
    );

    // Tests for Bls12_377
    batch_test!(
        [test_var_base_scalar_mul],
        [] => (
            Bls12_377,
            ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
