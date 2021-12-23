// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! Variable-base Curve Addition Gate

use crate::constraint_system::ecc::Point;
use crate::constraint_system::StandardComposer;
use ark_ec::models::twisted_edwards_extended::GroupAffine;
use ark_ec::models::TEModelParameters;
use ark_ec::PairingEngine;
use num_traits::{One, Zero};

impl<E, P> StandardComposer<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    /// Adds two curve points together using a curve addition gate
    /// Note that since the points are not fixed the generator is not a part of
    /// the circuit description, however it is less efficient for a program
    /// width of 4.
    pub fn point_addition_gate(
        &mut self,
        point_a: Point<E, P>,
        point_b: Point<E, P>,
    ) -> Point<E, P> {
        // In order to verify that two points were correctly added
        // without going over a degree 4 polynomial, we will need
        // x_1, y_1, x_2, y_2
        // x_3, y_3,      x_1 * y_2

        let x_1 = point_a.x;
        let y_1 = point_a.y;
        let x_2 = point_b.x;
        let y_2 = point_b.y;

        // Compute the resulting point
        let x_1_scalar = self.variables.get(&x_1).unwrap();
        let y_1_scalar = self.variables.get(&y_1).unwrap();
        let x_2_scalar = self.variables.get(&x_2).unwrap();
        let y_2_scalar = self.variables.get(&y_2).unwrap();

        let p1 = GroupAffine::<P>::new(*x_1_scalar, *y_1_scalar);
        let p2 = GroupAffine::<P>::new(*x_2_scalar, *y_2_scalar);

        let point = p1 + p2;
        let x_3_scalar = point.x;
        let y_3_scalar = point.y;

        let x1_scalar_y2_scalar = *x_1_scalar * y_2_scalar;

        // Add the rest of the prepared points into the composer
        let x_1_y_2 = self.add_input(x1_scalar_y2_scalar);
        let x_3 = self.add_input(x_3_scalar);
        let y_3 = self.add_input(y_3_scalar);

        self.w_l.extend(&[x_1, x_3]);
        self.w_r.extend(&[y_1, y_3]);
        self.w_o.extend(&[x_2, self.zero_var]);
        self.w_4.extend(&[y_2, x_1_y_2]);
        let zeros = [E::Fr::zero(), E::Fr::zero()];

        self.q_l.extend(&zeros);
        self.q_r.extend(&zeros);
        self.q_c.extend(&zeros);
        self.q_o.extend(&zeros);
        self.q_m.extend(&zeros);
        self.q_4.extend(&zeros);
        self.q_arith.extend(&zeros);
        self.q_range.extend(&zeros);
        self.q_logic.extend(&zeros);
        self.q_fixed_group_add.extend(&zeros);

        self.q_variable_group_add.push(E::Fr::one());
        self.q_variable_group_add.push(E::Fr::zero());

        self.perm.add_variables_to_map(x_1, y_1, x_2, y_2, self.n);
        self.n += 1;

        self.perm.add_variables_to_map(
            x_3,
            y_3,
            self.zero_var,
            x_1_y_2,
            self.n,
        );
        self.n += 1;

        Point::<E, P>::new(x_3, y_3)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{batch_test, constraint_system::helper::*};
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_ff::Field;

    /// Adds two curve points together using the classical point addition
    /// algorithm. This method is slower than WNAF and is just meant to be the
    /// source of truth to test the WNAF method.
    pub fn classical_point_addition<E, P>(
        composer: &mut StandardComposer<E, P>,
        point_a: Point<E, P>,
        point_b: Point<E, P>,
    ) -> Point<E, P>
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let zero = composer.zero_var;
        let x1 = point_a.x;
        let y1 = point_a.y;

        let x2 = point_b.x;
        let y2 = point_b.y;

        // x1 * y2
        let x1_y2 = composer.arithmetic_gate(|gate| {
            gate.mul(E::Fr::one()).witness(x1, y2, None)
        });
        // y1 * x2
        let y1_x2 = composer.arithmetic_gate(|gate| {
            gate.mul(E::Fr::one()).witness(y1, x2, None)
        });
        // y1 * y2
        let y1_y2 = composer.arithmetic_gate(|gate| {
            gate.mul(E::Fr::one()).witness(y1, y2, None)
        });
        // x1 * x2
        let x1_x2 = composer.arithmetic_gate(|gate| {
            gate.mul(E::Fr::one()).witness(x1, x2, None)
        });
        // d x1x2 * y1y2
        let d_x1_x2_y1_y2 = composer.arithmetic_gate(|gate| {
            gate.mul(P::COEFF_D).witness(x1_x2, y1_y2, None)
        });

        // x1y2 + y1x2
        let x_numerator = composer.arithmetic_gate(|gate| {
            gate.witness(x1_y2, y1_x2, None)
                .add(E::Fr::one(), E::Fr::one())
        });

        // y1y2 - a * x1x2 (a=-1) => y1y2 + x1x2
        let y_numerator = composer.arithmetic_gate(|gate| {
            gate.witness(y1_y2, x1_x2, None)
                .add(E::Fr::one(), E::Fr::one())
        });

        // 1 + dx1x2y1y2
        let x_denominator = composer.arithmetic_gate(|gate| {
            gate.witness(d_x1_x2_y1_y2, zero, None)
                .add(E::Fr::one(), E::Fr::zero())
                .constant(E::Fr::one())
        });

        // Compute the inverse
        let inv_x_denom = composer
            .variables
            .get(&x_denominator)
            .unwrap()
            .inverse()
            .unwrap();
        let inv_x_denom = composer.add_input(inv_x_denom);

        // Assert that we actually have the inverse
        // inv_x * x = 1
        composer.arithmetic_gate(|gate| {
            gate.witness(x_denominator, inv_x_denom, Some(zero))
                .mul(E::Fr::one())
                .constant(-E::Fr::one())
        });

        // 1 - dx1x2y1y2
        let y_denominator = composer.arithmetic_gate(|gate| {
            gate.witness(d_x1_x2_y1_y2, zero, None)
                .add(-E::Fr::one(), E::Fr::zero())
                .constant(E::Fr::one())
        });

        let inv_y_denom = composer
            .variables
            .get(&y_denominator)
            .unwrap()
            .inverse()
            .unwrap();

        let inv_y_denom = composer.add_input(inv_y_denom);
        // Assert that we actually have the inverse
        // inv_y * y = 1
        composer.arithmetic_gate(|gate| {
            gate.mul(E::Fr::one())
                .witness(y_denominator, inv_y_denom, Some(zero))
                .constant(-E::Fr::one())
        });

        // We can now use the inverses

        let x_3 = composer.arithmetic_gate(|gate| {
            gate.mul(E::Fr::one())
                .witness(inv_x_denom, x_numerator, None)
        });

        let y_3 = composer.arithmetic_gate(|gate| {
            gate.mul(E::Fr::one())
                .witness(inv_y_denom, y_numerator, None)
        });

        Point::new(x_3, y_3)
    }

    fn test_curve_addition<E, P>()
    where
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    {
        let res = gadget_tester(
            |composer: &mut StandardComposer<E, P>| {
                let (x, y) = P::AFFINE_GENERATOR_COEFFS;
                let generator = GroupAffine::new(x, y);
                let x_var = composer.add_input(x);
                let y_var = composer.add_input(y);
                let expected_point = generator + generator;
                let point_a = Point::new(x_var, y_var);
                let point_b = Point::new(x_var, y_var);

                let point = composer.point_addition_gate(point_a, point_b);
                let point2 =
                    classical_point_addition(composer, point_a, point_b);

                composer.assert_equal_point(point, point2);

                composer.assert_equal_public_point(point, expected_point);
            },
            2000,
        );
        assert!(res.is_ok());
    }

    batch_test!(
        [test_curve_addition],
        []
        => (
            Bls12_381,
            ark_ed_on_bls12_381::EdwardsParameters
        )
    );

    batch_test!(
        [test_curve_addition],
        [] => (
            Bls12_377,
            ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
