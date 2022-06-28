// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Variable-base Curve Addition Gate

use crate::{
    constraint_system::{ecc::Point, StandardComposer},
    parameters::CircuitParameters,
    proof_system::ecc::SWEmbeddedCurve,
};
use ark_ec::models::{
    short_weierstrass_jacobian::GroupAffine as SWGroupAffine, SWModelParameters,
};
use ark_ff::PrimeField;

impl<P, EmbeddedBaseField, EmbeddedCurveParameters> StandardComposer<P>
where
    EmbeddedBaseField: PrimeField,
    EmbeddedCurveParameters: SWModelParameters<BaseField = EmbeddedBaseField>,
    P: CircuitParameters<
        ScalarField = EmbeddedBaseField,
        EmbeddedCurve = SWEmbeddedCurve<P>,
        EmbeddedCurveParameters = EmbeddedCurveParameters,
    >,
{
    /// Adds two curve points together using a curve addition gate
    /// Note that since the points are not fixed the generator is not a part of
    /// the circuit description, however it is less efficient for a program
    /// width of 4.
    pub fn point_addition_gate_sw(
        &mut self,
        point_a: Point,
        point_b: Point,
    ) -> Point {
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

        let p1 = SWGroupAffine::<EmbeddedCurveParameters>::new(
            *x_1_scalar,
            *y_1_scalar,
            false,
        );
        let p2 = SWGroupAffine::<EmbeddedCurveParameters>::new(
            *x_2_scalar,
            *y_2_scalar,
            false,
        );

        let point = p1 + p2;
        let x_3_scalar = point.x;
        let y_3_scalar = point.y;

        // Add the rest of the prepared points into the composer
        let x_3 = self.add_input(x_3_scalar);
        let y_3 = self.add_input(y_3_scalar);

        self.w_l.extend(&[x_1, x_3]);
        self.w_r.extend(&[y_1, y_3]);
        self.w_o.extend(&[x_2, self.zero_var]);
        self.w_4.extend(&[y_2, self.zero_var]);
        let zeros = [P::ScalarField::zero(), P::ScalarField::zero()];

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
        self.q_lookup.extend(&zeros);

        self.q_variable_group_add.push(P::ScalarField::one());
        self.q_variable_group_add.push(P::ScalarField::zero());

        self.perm.add_variables_to_map(x_1, y_1, x_2, y_2, self.n);
        self.n += 1;

        self.perm.add_variables_to_map(
            x_3,
            y_3,
            self.zero_var,
            self.zero_var,
            self.n,
        );
        self.n += 1;

        Point::new(x_3, y_3)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        batch_test_embedded,
        constraint_system::helper::*,
        parameters::{test::*, CircuitParameters},
        proof_system::ecc::SWEmbeddedCurve,
    };

    fn test_curve_addition<
        CircuitParams,
        EmbeddedBaseField,
        EmbeddedCurveParameters,
    >()
    where
        EmbeddedBaseField: PrimeField,
        EmbeddedCurveParameters:
            SWModelParameters<BaseField = EmbeddedBaseField>,
        CircuitParams: CircuitParameters<
            ScalarField = EmbeddedBaseField,
            EmbeddedCurve = SWEmbeddedCurve<CircuitParams>,
            EmbeddedCurveParameters = EmbeddedCurveParameters,
        >,
    {
        let res = gadget_tester::<CircuitParams>(
            |composer: &mut StandardComposer<CircuitParams>| {
                let (x, y) = EmbeddedCurveParameters::AFFINE_GENERATOR_COEFFS;
                let generator =
                    SWGroupAffine::<EmbeddedCurveParameters>::new(x, y, false);
                let x_var = composer.add_input(x);
                let y_var = composer.add_input(y);
                let expected_point = generator + generator;
                let point_a = Point::new(x_var, y_var);
                let point_b = Point::new(x_var, y_var);

                let point = composer.point_addition_gate_sw(point_a, point_b);

                composer.assert_equal_public_point_sw(point, expected_point);
            },
            2000,
        );
        assert!(res.is_ok());
    }

    batch_test_embedded!(
        [test_curve_addition],
        []
        => [ BW6_761_KZG, Pallas_IPA, Vesta_IPA ]
    );
}
