// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! Range Gate

use crate::proof_system::GateConstraint;
use crate::proof_system::GateValues;
use ark_ff::Field;
use core::marker::PhantomData;

/// Range Gate
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct Range<F>(PhantomData<F>)
where
    F: Field;

impl<F> GateConstraint<F> for Range<F>
where
    F: Field,
{
    #[inline]
    fn constraints(separation_challenge: F, values: GateValues<F>) -> F {
        let four = F::from(4u64);
        let kappa = separation_challenge.square();
        let kappa_sq = kappa.square();
        let kappa_cu = kappa_sq * kappa;
        let b_1 = delta(values.output - four * values.fourth);
        let b_2 = delta(values.right - four * values.output) * kappa;
        let b_3 = delta(values.left - four * values.right) * kappa_sq;
        let b_4 = delta(values.fourth_next - four * values.left) * kappa_cu;
        (b_1 + b_2 + b_3 + b_4) * separation_challenge
    }
}

/// Computes `f(f-1)(f-2)(f-3)`.
fn delta<F>(f: F) -> F
where
    F: Field,
{
    let f_1 = f - F::one();
    let f_2 = f - F::from(2_u64);
    let f_3 = f - F::from(3_u64);
    f * f_1 * f_2 * f_3
}
