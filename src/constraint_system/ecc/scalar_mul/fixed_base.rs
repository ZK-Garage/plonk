// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::constraint_system::ecc::curve_addition::fixed_base_gate::WnafRound;
use crate::constraint_system::ecc::Point;
use crate::constraint_system::{variable::Variable, StandardComposer};
use ark_ec::models::twisted_edwards_extended::{GroupAffine, GroupProjective};
use ark_ec::models::TEModelParameters;
use ark_ec::{PairingEngine, ProjectiveCurve};
use ark_ff::{BigInteger, PrimeField};
use num_traits::{One, Zero};

fn compute_wnaf_point_multiples<P: TEModelParameters>() -> Vec<GroupAffine<P>> {
    let mut multiples = vec![
        GroupProjective::<P>::default();
        <P::BaseField as PrimeField>::Params::MODULUS_BITS
            as usize
    ];
    let (x, y) = P::AFFINE_GENERATOR_COEFFS;
    multiples[0] = GroupAffine::new(x, y).into();
    for i in 1..<P::BaseField as PrimeField>::Params::MODULUS_BITS as usize {
        multiples[i] = multiples[i - 1].double();
    }

    ProjectiveCurve::batch_normalization_into_affine(&mut multiples)
}

impl<
        E: PairingEngine<Fr = T::BaseField>,
        T: ProjectiveCurve<BaseField = P::BaseField>,
        P: TEModelParameters,
    > StandardComposer<E, T, P>
{
    /// Adds an elliptic curve Scalar multiplication gate to the circuit
    /// description.
    ///
    /// # Note
    /// This function is optimized for fixed base ops **ONLY** and therefore,
    /// the **ONLY** `generator` inputs that should be passed to this
    /// function as inputs are [`dusk_jubjub::GENERATOR`] or
    /// [`dusk_jubjub::GENERATOR_NUMS`].
    pub fn fixed_base_scalar_mul(
        &mut self,
        jubjub_scalar: Variable,
    ) -> Point<E, T, P> {
        let num_bits =
            <P::BaseField as PrimeField>::Params::MODULUS_BITS as usize;
        // compute 2^iG
        let mut point_multiples = compute_wnaf_point_multiples();
        point_multiples.reverse();

        // Fetch the raw scalar value as bls scalar, then convert to a jubjub
        // scalar
        // XXX: Not very Tidy, impl From function in JubJub
        let raw_bls_scalar = self.variables.get(&jubjub_scalar).unwrap();
        let raw_jubjub_scalar =
            <P::BaseField as PrimeField>::from_le_bytes_mod_order(
                &raw_bls_scalar.to_bytes_le(),
            )
            .unwrap();

        // Convert scalar to wnaf_2(k)
        let wnaf_entries = raw_jubjub_scalar.compute_windowed_naf(2);
        assert_eq!(wnaf_entries.len(), num_bits);

        // Initialise the accumulators
        let mut scalar_acc = vec![E::Fr::zero()];
        let mut point_acc = vec![GroupAffine::<P>::identity()];

        // Auxillary point to help with checks on the backend
        let mut xy_alphas = Vec::new();

        // Load values into accumulators based on wnaf entries
        for (i, entry) in wnaf_entries.iter().rev().enumerate() {
            // Based on the WNAF, we decide what scalar and point to add
            let (scalar_to_add, point_to_add) = match entry {
            0 => { (E::Fr::zero(), GroupAffine::<P>::identity())},
            -1 => {(E::Fr::one().neg(), -point_multiples[i])},
            1 => {(E::Fr::one(), point_multiples[i])},
            _ => unreachable!("Currently WNAF_2(k) is supported. The possible values are 1, -1 and 0. Current entry is {}", entry),
        };

            let prev_accumulator = E::Fr::from(2u64) * scalar_acc[i];
            scalar_acc.push(prev_accumulator + scalar_to_add);
            point_acc.push(
                (GroupProjective::<P>::from(point_acc[i])
                    + GroupProjective::<P>::from(point_to_add))
                .into(),
            );

            let x_alpha = point_to_add.get_x();
            let y_alpha = point_to_add.get_y();

            xy_alphas.push(x_alpha * y_alpha);
        }

        for i in 0..num_bits {
            let acc_x = self.add_input(point_acc[i].get_x());
            let acc_y = self.add_input(point_acc[i].get_y());

            let accumulated_bit = self.add_input(scalar_acc[i]);

            // We constrain the point accumulator to start from the Identity
            // point and the Scalar accumulator to start from zero
            if i == 0 {
                self.constrain_to_constant(acc_x, E::Fr::zero(), None);
                self.constrain_to_constant(acc_y, E::Fr::one(), None);
                self.constrain_to_constant(
                    accumulated_bit,
                    E::Fr::zero(),
                    None,
                );
            }

            let x_beta = point_multiples[i].get_x();
            let y_beta = point_multiples[i].get_y();

            let xy_alpha = self.add_input(xy_alphas[i]);

            let xy_beta = x_beta * y_beta;

            let wnaf_round = WnafRound::new_wnaf(
                acc_x,
                acc_y,
                accumulated_bit,
                xy_alpha,
                x_beta,
                y_beta,
                xy_beta,
            );

            self.fixed_group_add(wnaf_round);
        }

        // Add last gate, but do not activate it for ECC
        // It is for use with the previous gate
        let acc_x = self.add_input(point_acc[num_bits].get_x());
        let acc_y = self.add_input(point_acc[num_bits].get_y());
        let xy_alpha = self.zero_var;
        let last_accumulated_bit = self.add_input(scalar_acc[num_bits]);

        self.big_add_gate(
            acc_x,
            acc_y,
            xy_alpha,
            Some(last_accumulated_bit),
            E::Fr::zero(),
            E::Fr::zero(),
            E::Fr::zero(),
            E::Fr::zero(),
            E::Fr::zero(),
            None,
        );

        // Constrain the last element in the accumulator to be equal to the
        // input jubjub scalar
        self.assert_equal(last_accumulated_bit, jubjub_scalar);

        Point::<E, T, P>::new(acc_x, acc_y)
    }
}

/*
#[cfg(feature = "std")]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::constraint_system::helper::*;
    use ark_bls12_381::Fr as BlsScalar;
    use ark_ed_on_bls12_381::{EdwardsAffine, EdwardsParameters, Fr};

    #[test]
    fn test_ecc_constraint() {
        let res = gadget_tester(
            |composer| {
                let scalar = Fr::from_le_bytes_mod_order(&[
                    182, 44, 247, 214, 94, 14, 151, 208, 130, 16, 200, 204,
                    147, 32, 104, 166, 0, 59, 52, 1, 1, 59, 103, 6, 169, 175,
                    51, 101, 234, 180, 125, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0,
                ]);
                let bls_scalar =
                    BlsScalar::from_bytes(&scalar.to_bytes()).unwrap();
                let secret_scalar = composer.add_input(bls_scalar);

                let (x, y) = EdwardsParameters::AFFINE_GENERATOR_COEFFS;
                let generator = ark_ed_on_bls12_381::EdwardsAffine::new(x, y);
                let expected_point: EdwardsAffine = (generator * scalar).into();

                let point_scalar =
                    composer.fixed_base_scalar_mul(secret_scalar);

                composer
                    .assert_equal_public_point(point_scalar, expected_point);
            },
            600,
        );
        assert!(res.is_ok());
    }

    #[test]
    fn test_ecc_constraint_zero() {
        let res = gadget_tester(
            |composer| {
                let scalar = Fr::zero();
                let bls_scalar =
                    BlsScalar::from_le_bytes_mod_order(&scalar.to_bytes())
                        .unwrap();
                let secret_scalar = composer.add_input(bls_scalar);

                let (x, y) = EdwardsParameters::AFFINE_GENERATOR_COEFFS;
                let generator = ark_ed_on_bls12_381::EdwardsAffine::new(x, y);
                let expected_point: EdwardsAffine = (generator * scalar).into();

                let point_scalar =
                    composer.fixed_base_scalar_mul(secret_scalar);

                composer
                    .assert_equal_public_point(point_scalar, expected_point);
            },
            600,
        );
        assert!(res.is_ok());
    }
    #[test]
    fn test_ecc_constraint_should_fail() {
        let res = gadget_tester(
            |composer| {
                let scalar = Fr::from(100u64);
                let bls_scalar =
                    BlsScalar::from_le_bytes_mod_order(&scalar.to_bytes())
                        .unwrap();
                let secret_scalar = composer.add_input(bls_scalar);
                // Fails because we are not multiplying by the GENERATOR, it is
                // double
                let (x, y) = EdwardsParameters::AFFINE_GENERATOR_COEFFS;
                let generator = ark_ed_on_bls12_381::EdwardsAffine::new(x, y);
                let double_gen = generator.double();

                let expected_point: GroupAffine<EdwardsParameters> =
                    (double_gen * scalar).into();

                let point_scalar =
                    composer.fixed_base_scalar_mul(secret_scalar);

                composer
                    .assert_equal_public_point(point_scalar, expected_point);
            },
            600,
        );

        assert!(res.is_err());
    }
    #[test]
    fn test_point_addition() {
        let res = gadget_tester(
            |composer| {
                let (x, y) = EdwardsParameters::AFFINE_GENERATOR_COEFFS;
                let generator = ark_ed_on_bls12_381::EdwardsAffine::new(x, y);

                let point_a = generator;
                let point_b = point_a.double();
                let expected_point = point_a + point_b;

                let affine_point_a: GroupAffine<EdwardsParameters> =
                    point_a.into();
                let affine_point_b: GroupAffine<EdwardsParameters> =
                    point_b.into();
                let affine_expected_point: GroupAffine<EdwardsParameters> =
                    expected_point.into();

                let var_point_a_x = composer.add_input(affine_point_a.get_x());
                let var_point_a_y = composer.add_input(affine_point_a.get_y());
                let point_a =
                    Point::<E, T, P>::new(var_point_a_x, var_point_a_y);
                let var_point_b_x = composer.add_input(affine_point_b.get_x());
                let var_point_b_y = composer.add_input(affine_point_b.get_y());
                let point_b =
                    Point::<E, T, P>::new(var_point_b_x, var_point_b_y);
                let new_point = composer.point_addition_gate(point_a, point_b);

                composer.assert_equal_public_point(
                    new_point,
                    affine_expected_point,
                );
            },
            600,
        );

        assert!(res.is_ok());
    }
    #[test]
    #[allow(non_snake_case)]
    fn test_pedersen_hash() {
        let res = gadget_tester(
            |composer| {
                let (x, y) = EdwardsParameters::AFFINE_GENERATOR_COEFFS;
                let generator = ark_ed_on_bls12_381::EdwardsAffine::new(x, y);
                // First component
                let scalar_a = Fr::from(112233u64);
                let bls_scalar =
                    BlsScalar::from_le_bytes_mod_order(&scalar_a.to_bytes())
                        .unwrap();
                let secret_scalar_a = composer.add_input(bls_scalar);
                let point_a = generator;
                let c_a: GroupAffine<EdwardsParameters> =
                    (point_a * scalar_a).into();

                // Second component
                let scalar_b = Fr::from(445566u64);
                let bls_scalar =
                    BlsScalar::from_le_bytes_mod_order(&scalar_b.to_bytes())
                        .unwrap();
                let secret_scalar_b = composer.add_input(bls_scalar);
                let point_b = point_a.double() + point_a;
                let c_b: GroupAffine<EdwardsParameters> =
                    (point_b * scalar_b).into();

                // Expected pedersen hash
                let expected_point: GroupAffine<EdwardsParameters> =
                    (point_a * scalar_a + point_b * scalar_b).into();

                // To check this pedersen commitment, we will need to do:
                // - Two scalar multiplications
                // - One curve addition
                //
                // Scalar multiplications
                let aG = composer.fixed_base_scalar_mul(secret_scalar_a);
                let bH = composer.fixed_base_scalar_mul(secret_scalar_b);

                // Depending on the context, one can check if the resulting aG
                // and bH are as expected
                //
                composer.assert_equal_public_point(aG, c_a);
                composer.assert_equal_public_point(bH, c_b);

                // Curve addition
                let commitment = composer.point_addition_gate(aG, bH);

                // Add final constraints to ensure that the commitment that we
                // computed is equal to the public point
                composer.assert_equal_public_point(commitment, expected_point);
            },
            1024,
        );
        assert!(res.is_ok());
    }
    #[test]
    #[allow(non_snake_case)]
    fn test_pedersen_balance() {
        let res = gadget_tester(
            |composer| {
                // First component
                let scalar_a = Fr::from(25u64);
                let bls_scalar_a =
                    BlsScalar::from_le_bytes_mod_order(&scalar_a.to_bytes())
                        .unwrap();
                let secret_scalar_a = composer.add_input(bls_scalar_a);
                // Second component
                let scalar_b = Fr::from(30u64);
                let bls_scalar_b =
                    BlsScalar::from_le_bytes_mod_order(&scalar_b.to_bytes())
                        .unwrap();
                let secret_scalar_b = composer.add_input(bls_scalar_b);
                // Third component
                let scalar_c = Fr::from(10u64);
                let bls_scalar_c =
                    BlsScalar::from_le_bytes_mod_order(&scalar_c.to_bytes())
                        .unwrap();
                let secret_scalar_c = composer.add_input(bls_scalar_c);
                // Fourth component
                let scalar_d = Fr::from(45u64);
                let bls_scalar_d =
                    BlsScalar::from_le_bytes_mod_order(&scalar_d.to_bytes())
                        .unwrap();
                let secret_scalar_d = composer.add_input(bls_scalar_d);

                let (x, y) = EdwardsParameters::AFFINE_GENERATOR_COEFFS;
                let gen = ark_ed_on_bls12_381::EdwardsAffine::new(x, y);

                let expected_lhs: GroupAffine<EdwardsParameters> =
                    (gen * (scalar_a + scalar_b)).into();
                let expected_rhs: GroupAffine<EdwardsParameters> =
                    (gen * (scalar_c + scalar_d)).into();

                let P1 = composer.fixed_base_scalar_mul(secret_scalar_a);
                let P2 = composer.fixed_base_scalar_mul(secret_scalar_b);
                let P3 = composer.fixed_base_scalar_mul(secret_scalar_c);
                let P4 = composer.fixed_base_scalar_mul(secret_scalar_d);

                let commitment_a = composer.point_addition_gate(P1, P2);
                let commitment_b = composer.point_addition_gate(P3, P4);

                composer.assert_equal_point(commitment_a, commitment_b);

                composer.assert_equal_public_point(commitment_a, expected_lhs);
                composer.assert_equal_public_point(commitment_b, expected_rhs);
            },
            2048,
        );
        assert!(res.is_ok());
    }
}

*/
