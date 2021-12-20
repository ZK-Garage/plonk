// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-INFRA. All rights reserved.

//! Fixed-Base Curve Addition Gate

use crate::constraint_system::StandardComposer;
use crate::constraint_system::Variable;
use ark_ec::models::TEModelParameters;
use ark_ec::PairingEngine;
use core::marker::PhantomData;
use num_traits::{One, Zero};

/// Contains all of the components needed to verify that a bit scalar
/// multiplication was computed correctly.
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug)]
pub struct WnafRound<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    /// This is the accumulated x coordinate point that we wish to add (so
    /// far, it depends on where you are in the scalar mul) it is linked to
    /// the wnaf entry, so must not be revealed
    pub acc_x: Variable,

    /// This is the accumulated y coordinate
    pub acc_y: Variable,

    /// This is the wnaf accumulated entry
    /// For all intents and purposes, you can think of this as the secret bit
    pub accumulated_bit: Variable,

    /// This is the multiplication of x_\alpha * y_\alpha
    /// we need this as a distinct wire, so that the degree of the polynomial
    /// does not go over 4
    pub xy_alpha: Variable,

    /// This is the possible x co-ordinate of the wnaf point we are going to
    /// add Actual x-co-ordinate = b_i * x_\beta
    pub x_beta: P::BaseField,

    /// This is the possible y co-ordinate of the wnaf point we are going to
    /// add Actual y coordinate = (b_i)^2 [y_\beta -1] + 1
    pub y_beta: P::BaseField,

    /// This is the multiplication of x_\beta * y_\beta
    pub xy_beta: P::BaseField,

    /// Type Parameter Marker
    __: PhantomData<E>,
}

impl<E, P> StandardComposer<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    /// Generates a new structure for preparing a [`WnafRound`] ROUND.
    pub(crate) fn new_wnaf(
        acc_x: Variable,
        acc_y: Variable,
        accumulated_bit: Variable,
        xy_alpha: Variable,
        x_beta: P::BaseField,
        y_beta: P::BaseField,
        xy_beta: P::BaseField,
    ) -> WnafRound<E, P> {
        WnafRound {
            acc_x,
            acc_y,
            accumulated_bit,
            xy_alpha,
            x_beta,
            y_beta,
            xy_beta,
            __: PhantomData,
        }
    }

    /// Fixed group addition of a point.
    pub(crate) fn fixed_group_add(&mut self, wnaf_round: WnafRound<E, P>) {
        self.w_l.push(wnaf_round.acc_x);
        self.w_r.push(wnaf_round.acc_y);
        self.w_o.push(wnaf_round.xy_alpha);
        self.w_4.push(wnaf_round.accumulated_bit);

        self.q_l.push(wnaf_round.x_beta);
        self.q_r.push(wnaf_round.y_beta);

        self.q_c.push(wnaf_round.xy_beta);
        self.q_o.push(E::Fr::zero());
        self.q_fixed_group_add.push(E::Fr::one());
        self.q_variable_group_add.push(E::Fr::zero());

        self.q_m.push(E::Fr::zero());
        self.q_4.push(E::Fr::zero());
        self.q_arith.push(E::Fr::zero());
        self.q_range.push(E::Fr::zero());
        self.q_logic.push(E::Fr::zero());

        self.perm.add_variables_to_map(
            wnaf_round.acc_x,
            wnaf_round.acc_y,
            wnaf_round.xy_alpha,
            wnaf_round.accumulated_bit,
            self.n,
        );

        self.n += 1;
    }
}
