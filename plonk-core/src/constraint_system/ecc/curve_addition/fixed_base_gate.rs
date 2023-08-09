// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Fixed-Base Curve Addition Gate

use crate::constraint_system::{StandardComposer, Variable};
use ark_ec::models::TEModelParameters;
use ark_ff::PrimeField;

/// Contains all of the components needed to verify that a bit scalar
/// multiplication was computed correctly.
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug)]
pub struct WnafRound<P>
where
    P: TEModelParameters,
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
}

impl<F, P> StandardComposer<F, P>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    /// Generates a new structure for preparing a [`WnafRound`] ROUND.
    pub(crate) fn new_wnaf(
        acc_x: Variable,
        acc_y: Variable,
        accumulated_bit: Variable,
        xy_alpha: Variable,
        x_beta: F,
        y_beta: F,
        xy_beta: F,
    ) -> WnafRound<P> {
        WnafRound {
            acc_x,
            acc_y,
            accumulated_bit,
            xy_alpha,
            x_beta,
            y_beta,
            xy_beta,
        }
    }

    /// Fixed group addition of a point.
    pub(crate) fn fixed_group_add(&mut self, wnaf_round: WnafRound<P>) {
        self.w_l.push(wnaf_round.acc_x);
        self.w_r.push(wnaf_round.acc_y);
        self.w_o.push(wnaf_round.xy_alpha);
        self.w_4.push(wnaf_round.accumulated_bit);

        self.q_l.push(wnaf_round.x_beta);
        self.q_r.push(wnaf_round.y_beta);

        self.q_c.push(wnaf_round.xy_beta);
        self.q_o.push(F::zero());
        self.q_fixed_group_add.push(F::one());
        self.q_variable_group_add.push(F::zero());

        self.q_m.push(F::zero());
        self.q_4.push(F::zero());
        self.q_arith.push(F::zero());
        self.q_range.push(F::zero());
        self.q_logic.push(F::zero());
        self.q_lookup.push(F::zero());

        // add high degree selectors
        self.q_hl.push(F::zero());
        self.q_hr.push(F::zero());
        self.q_h4.push(F::zero());

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
