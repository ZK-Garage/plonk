// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

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
    /// bit == 1 => lhs,
    /// bit == 0 => rhs,
    /// ```
    ///
    /// # Note
    ///
    /// The `bit` used as input which is a [`Variable`] should have previously
    /// been constrained to be either `1` or `0` using a boolean constraint.
    /// See: [`StandardComposer::boolean_gate`].
    pub fn conditional_point_select(
        &mut self,
        lhs: Point<E, P>,
        rhs: Point<E, P>,
        bit: Variable,
    ) -> Point<E, P> {
        Point::new(
            self.conditional_select(bit, lhs.x, rhs.x),
            self.conditional_select(bit, lhs.y, rhs.y),
        )
    }

    /// Adds to the circuit description the conditional selection of the
    /// identity point:
    ///
    /// ```text
    /// bit == 1 => value,
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

/* TODO: Do we need this test?
#[cfg(feature = "std")]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraint_system::helper::*;
    use ark_bls12_381::Fr as BlsScalar;

    #[test]
    fn test_conditional_select_point() {
        let res = gadget_tester(
            |composer| {
                let bit_1 = composer.add_input(BlsScalar::one());
                let bit_0 = composer.zero_var();

                let point_a = Point::identity(composer);
                let point_b = Point {
                    x: composer.add_input(BlsScalar::from(10u64)),
                    y: composer.add_input(BlsScalar::from(20u64)),
                };

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
}
*/
