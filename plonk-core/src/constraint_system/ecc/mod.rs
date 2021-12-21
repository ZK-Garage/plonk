// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! Elliptic Curve Gates

pub mod curve_addition;
pub mod scalar_mul;

use crate::constraint_system::{variable::Variable, StandardComposer};
use ark_ec::{
    twisted_edwards_extended::GroupAffine, PairingEngine, TEModelParameters,
};
use core::marker::PhantomData;
use num_traits::{One, Zero};

/// Represents a point of the embeded curve in the circuit
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug)]
pub struct Point<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    /// `X`-coordinate
    x: Variable,

    /// `Y`-coordinate
    y: Variable,

    /// Type Parameter Marker
    __: PhantomData<(E, P)>,
}

impl<E, P> Point<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    /// Builds a new [`Point`] from `X` and `Y` coordinates.
    ///
    /// # Safety
    ///
    /// This method should only be called when we have a guarantee in some
    /// [`StandardComposer`] that these two are valid coordinates for an
    /// elliptic curve point in affine form.
    pub fn new(x: Variable, y: Variable) -> Self {
        Self {
            x,
            y,
            __: PhantomData,
        }
    }

    /// Returns an identity point.
    pub fn identity(composer: &mut StandardComposer<E, P>) -> Self {
        let one = composer.add_witness_to_circuit_description(E::Fr::one());
        Self::new(composer.zero_var, one)
    }

    /// Returns the `X`-coordinate of `self`.
    pub fn x(&self) -> &Variable {
        &self.x
    }

    /// Returns the `Y`-coordinate of `self`.
    pub fn y(&self) -> &Variable {
        &self.y
    }
}

impl<E, P> StandardComposer<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    /// Converts an embeded curve point into a constraint system Point
    /// without constraining the values
    pub fn add_affine(&mut self, affine: GroupAffine<P>) -> Point<E, P> {
        Point::new(self.add_input(affine.x), self.add_input(affine.y))
    }

    /// Converts an embeded curve point into a constraint system Point
    /// without constraining the values
    pub fn add_public_affine(&mut self, affine: GroupAffine<P>) -> Point<E, P> {
        let point = self.add_affine(affine);
        self.constrain_to_constant(point.x, E::Fr::zero(), Some(-affine.x));
        self.constrain_to_constant(point.y, E::Fr::zero(), Some(-affine.y));
        point
    }

    /// Add the provided affine point as a circuit description and return its
    /// constrained witness value
    pub fn add_affine_to_circuit_description(
        &mut self,
        affine: GroupAffine<P>,
    ) -> Point<E, P> {
        // NOTE: Not using individual gates because one of these may be zero.
        Point::new(
            self.add_witness_to_circuit_description(affine.x),
            self.add_witness_to_circuit_description(affine.y),
        )
    }

    /// Asserts that a [`Point`] in the circuit is equal to a known public
    /// point.
    pub fn assert_equal_public_point(
        &mut self,
        point: Point<E, P>,
        public_point: GroupAffine<P>,
    ) {
        self.constrain_to_constant(
            point.x,
            E::Fr::zero(),
            Some(-public_point.x),
        );
        self.constrain_to_constant(
            point.y,
            E::Fr::zero(),
            Some(-public_point.y),
        );
    }

    /// Asserts that a point in the circuit is equal to another point in the
    /// circuit.
    pub fn assert_equal_point(&mut self, lhs: Point<E, P>, rhs: Point<E, P>) {
        self.assert_equal(lhs.x, rhs.x);
        self.assert_equal(lhs.y, rhs.y);
    }

    /// Adds to the circuit description the conditional selection of the
    /// a point between two of them:
    ///
    /// ```text
    /// bit == 1 => point_1,
    /// bit == 0 => point_0,
    /// ```
    ///
    /// # Note
    ///
    /// The `bit` used as input which is a [`Variable`] should have previously
    /// been constrained to be either `1` or `0` using a boolean constraint.
    /// See: [`StandardComposer::boolean_gate`].
    pub fn conditional_point_select(
        &mut self,
        point_1: Point<E, P>,
        point_0: Point<E, P>,
        bit: Variable,
    ) -> Point<E, P> {
        Point::new(
            self.conditional_select(bit, point_1.x, point_0.x),
            self.conditional_select(bit, point_1.y, point_0.y),
        )
    }

    /// Adds to the circuit description the conditional negation of a point:
    /// bit == 1 => -value,
    /// bit == 0 => value,
    ///
    /// # Note
    /// The `bit` used as input which is a [`Variable`] should had previously
    /// been constrained to be either 1 or 0 using a bool constrain. See:
    /// [`StandardComposer::boolean_gate`].
    pub fn conditional_point_neg(
        &mut self,
        bit: Variable,
        point_b: Point<E, P>,
    ) -> Point<E, P> {
        let x = point_b.x;
        let y = point_b.y;

        // negation of point (x, y) is (-x, y)
        let x_neg = self.add(
            (-E::Fr::one(), x),
            (E::Fr::zero(), self.zero_var),
            E::Fr::zero(),
            None,
        );
        let x_updated = self.conditional_select(bit, x_neg, x);

        Point::new(x_updated, y)
    }

    /// Adds to the circuit description the conditional selection of the
    /// identity point:
    ///
    /// ```text
    /// bit == 1 => point,
    /// bit == 0 => 1,
    /// ```
    ///
    /// # Note
    ///
    /// The `bit` used as input which is a [`Variable`] should have previously
    /// been constrained to be either `1` or `0` using a boolean constraint.
    /// See: [`StandardComposer::boolean_gate`].
    fn conditional_select_identity(
        &mut self,
        bit: Variable,
        point: Point<E, P>,
    ) -> Point<E, P> {
        Point::new(
            self.conditional_select_zero(bit, point.x),
            self.conditional_select_one(bit, point.y),
        )
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{batch_test, constraint_system::helper::*};
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;

    fn test_conditional_select_point<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let bit_1 = composer.add_input(E::Fr::one());
                let bit_0 = composer.zero_var();

                let point_a = Point::identity(composer);
                let point_b = Point::new(
                    composer.add_input(E::Fr::from(10u64)),
                    composer.add_input(E::Fr::from(20u64)),
                );

                let choice =
                    composer.conditional_point_select(point_a, point_b, bit_1);

                composer.assert_equal_point(point_a, choice);

                let choice =
                    composer.conditional_point_select(point_a, point_b, bit_0);
                composer.assert_equal_point(point_b, choice);
            },
            32,
        );
        assert!(res.is_ok());
    }

    fn test_conditional_point_neg<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let bit_1 = composer.add_input(E::Fr::one());
                let bit_0 = composer.zero_var();

                let point =
                    GroupAffine::new(E::Fr::from(10u64), E::Fr::from(20u64));
                let point_var = Point::new(
                    composer.add_input(point.x),
                    composer.add_input(point.y),
                );

                let neg_point =
                    composer.conditional_point_neg(bit_1, point_var);
                composer.assert_equal_public_point(neg_point, -point);

                let non_neg_point =
                    composer.conditional_point_neg(bit_0, point_var);
                composer.assert_equal_public_point(non_neg_point, point);
            },
            32,
        )
        .expect("test failed");
    }

    // Bls12-381 tests
    batch_test!(
        [
            test_conditional_select_point,
            test_conditional_point_neg
        ],
        [] => (
            Bls12_381,
            ark_ed_on_bls12_381::EdwardsParameters
        )
    );

    // Bls12-377 tests
    batch_test!(
        [
            test_conditional_select_point,
            test_conditional_point_neg
        ],
        [] => (
            Bls12_377,
            ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
