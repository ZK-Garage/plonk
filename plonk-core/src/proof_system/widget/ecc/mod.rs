// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Elliptic Curve Cryptography Gates

mod curve_addition;
mod fixed_base_scalar_mul;

pub use curve_addition::*;
pub use fixed_base_scalar_mul::*;

use crate::{
    constraint_system::{ecc::Point, StandardComposer},
    parameters::{CircuitParameters, EmbeddedCurve},
};
use ark_ec::{SWModelParameters, TEModelParameters};
use ark_ff::PrimeField;
use core::marker::PhantomData;

/// An embedded curve in twisted edwards form
pub struct TEEmbeddedCurve<P>(PhantomData<P>)
where
    P: CircuitParameters<EmbeddedCurve = Self>;

impl<P, EmbeddedBaseField, EmbeddedCurveParameters> EmbeddedCurve<P>
    for TEEmbeddedCurve<P>
where
    EmbeddedBaseField: PrimeField,
    EmbeddedCurveParameters: TEModelParameters<BaseField = EmbeddedBaseField>,
    P: CircuitParameters<
        ScalarField = EmbeddedBaseField,
        EmbeddedCurve = TEEmbeddedCurve<P>,
        EmbeddedCurveParameters = EmbeddedCurveParameters,
    >,
{
    /// Returns an identity point.
    fn identity(composer: &mut StandardComposer<P>) -> Point {
        let one =
            composer.add_witness_to_circuit_description(P::ScalarField::one());
        Point::new(composer.zero_var, one)
    }
}

/// An embedded curve in short weierstrass form
pub struct SWEmbeddedCurve<P>(PhantomData<P>)
where
    P: CircuitParameters<EmbeddedCurve = Self>;

impl<P, EmbeddedBaseField, EmbeddedCurveParameters> EmbeddedCurve<P>
    for SWEmbeddedCurve<P>
where
    EmbeddedBaseField: PrimeField,
    EmbeddedCurveParameters: SWModelParameters<BaseField = EmbeddedBaseField>,
    P: CircuitParameters<
        ScalarField = EmbeddedBaseField,
        EmbeddedCurve = SWEmbeddedCurve<P>,
        EmbeddedCurveParameters = EmbeddedCurveParameters,
    >,
{
    /// Returns an identity point.
    fn identity(composer: &mut StandardComposer<P>) -> Point {
        Point::new(composer.zero_var, composer.zero_var)
    }
}
