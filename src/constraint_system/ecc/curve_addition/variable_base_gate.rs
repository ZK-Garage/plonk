// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::constraint_system::ecc::Point;
use crate::constraint_system::StandardComposer;
use ark_ec::models::twisted_edwards_extended::GroupAffine;
use ark_ec::models::TEModelParameters;
use ark_ec::{PairingEngine, ProjectiveCurve};
use num_traits::{One, Zero};

impl<
        E: PairingEngine,
        T: ProjectiveCurve<BaseField = E::Fr>,
        P: TEModelParameters<BaseField = E::Fr>,
    > StandardComposer<E, T, P>
{
    /// Adds two curve points together using a curve addition gate
    /// Note that since the points are not fixed the generator is not a part of
    /// the circuit description, however it is less efficient for a program
    /// width of 4.
    pub fn point_addition_gate(
        &mut self,
        point_a: Point<E, T, P>,
        point_b: Point<E, T, P>,
    ) -> Point<E, T, P> {
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

        Point::<E, T, P>::new(x_3, y_3)
    }
}

#[cfg(test)]
mod variable_base_gate_tests {
    use super::*;
    use crate::constraint_system::helper::*;
    use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
    use ark_ed_on_bls12_381::{
        EdwardsAffine as JubjubAffine, EdwardsParameters as JubjubParameters,
        EdwardsProjective as JubjubProjective,
    };
    use ark_ff::Field;
    use num_traits::{One, Zero};
    /// Adds two curve points together using the classical point addition
    /// algorithm. This method is slower than WNaf and is just meant to be the
    /// source of truth to test the WNaf method.
    pub fn classical_point_addition(
        composer: &mut StandardComposer<
            Bls12_381,
            JubjubProjective,
            JubjubParameters,
        >,
        point_a: Point<Bls12_381, JubjubProjective, JubjubParameters>,
        point_b: Point<Bls12_381, JubjubProjective, JubjubParameters>,
    ) -> Point<Bls12_381, JubjubProjective, JubjubParameters> {
        let x1 = point_a.x;
        let y1 = point_a.y;

        let x2 = point_b.x;
        let y2 = point_b.y;

        // x1 * y2
        let x1_y2 =
            composer.mul(BlsScalar::one(), x1, y2, BlsScalar::zero(), None);
        // y1 * x2
        let y1_x2 =
            composer.mul(BlsScalar::one(), y1, x2, BlsScalar::zero(), None);
        // y1 * y2
        let y1_y2 =
            composer.mul(BlsScalar::one(), y1, y2, BlsScalar::zero(), None);
        // x1 * x2
        let x1_x2 =
            composer.mul(BlsScalar::one(), x1, x2, BlsScalar::zero(), None);
        // d x1x2 * y1y2
        let d_x1_x2_y1_y2 = composer.mul(
            JubjubParameters::COEFF_D,
            x1_x2,
            y1_y2,
            BlsScalar::zero(),
            None,
        );

        // x1y2 + y1x2
        let x_numerator = composer.add(
            (BlsScalar::one(), x1_y2),
            (BlsScalar::one(), y1_x2),
            BlsScalar::zero(),
            None,
        );

        // y1y2 - a * x1x2 (a=-1) => y1y2 + x1x2
        let y_numerator = composer.add(
            (BlsScalar::one(), y1_y2),
            (BlsScalar::one(), x1_x2),
            BlsScalar::zero(),
            None,
        );

        // 1 + dx1x2y1y2
        let x_denominator = composer.add(
            (BlsScalar::one(), d_x1_x2_y1_y2),
            (BlsScalar::zero(), composer.zero_var),
            BlsScalar::one(),
            None,
        );

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
        composer.mul_gate(
            x_denominator,
            inv_x_denom,
            composer.zero_var,
            BlsScalar::one(),
            BlsScalar::zero(),
            -BlsScalar::one(),
            None,
        );

        // 1 - dx1x2y1y2
        let y_denominator = composer.add(
            (-BlsScalar::one(), d_x1_x2_y1_y2),
            (BlsScalar::zero(), composer.zero_var),
            BlsScalar::one(),
            None,
        );
        let inv_y_denom = composer
            .variables
            .get(&y_denominator)
            .unwrap()
            .inverse()
            .unwrap();
        let inv_y_denom = composer.add_input(inv_y_denom);
        // Assert that we actually have the inverse
        // inv_y * y = 1
        composer.mul_gate(
            y_denominator,
            inv_y_denom,
            composer.zero_var,
            BlsScalar::one(),
            BlsScalar::zero(),
            -BlsScalar::one(),
            None,
        );

        // We can now use the inverses

        let x_3 = composer.mul(
            BlsScalar::one(),
            inv_x_denom,
            x_numerator,
            BlsScalar::zero(),
            None,
        );
        let y_3 = composer.mul(
            BlsScalar::one(),
            inv_y_denom,
            y_numerator,
            BlsScalar::zero(),
            None,
        );

        Point::new(x_3, y_3)
    }

    #[test]
    fn test_curve_addition() {
        let res = gadget_tester(
            |composer| {
                let (x, y) = JubjubParameters::AFFINE_GENERATOR_COEFFS;
                let generator = JubjubAffine::new(x, y);
                let x_var = composer.add_input(x);
                let y_var = composer.add_input(y);
                let expected_point: JubjubAffine =
                    (generator + generator).into();
                let point_a = Point::new(x_var, y_var);
                let point_b = Point::new(x_var, y_var);

                let point = composer.point_addition_gate(point_a, point_b);
                let point2 =
                    classical_point_addition(composer, point_a, point_b);

                composer.assert_equal_point(point, point2);

                composer
                    .assert_equal_public_point(point.into(), expected_point);
            },
            2000,
        );
        assert!(res.is_ok());
    }
}
