// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use super::ProverKey;
use crate::proof_system::ecc::CurveAddition;
use crate::proof_system::ecc::FixedBaseScalarMul;
use crate::proof_system::logic::Logic;
use crate::proof_system::range::Range;
use crate::proof_system::widget::GateConstraint;
use crate::proof_system::GateValues;
use crate::util::EvaluationDomainExt;
use ark_ec::{PairingEngine, TEModelParameters};
use ark_ff::Field;
use ark_ff::PrimeField;
use ark_poly::{
    univariate::DensePolynomial, GeneralEvaluationDomain, Polynomial,
};
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write,
};

/// Evaluations at points `z` or and `z * root of unity`
#[allow(dead_code)]
pub(crate) struct Evaluations<F: PrimeField> {
    pub(crate) proof: ProofEvaluations<F>,
    // Evaluation of the linearisation sigma polynomial at `z`
    pub(crate) quot_eval: F,
}

/// Subset of all of the evaluations. These evaluations
/// are added to the [`Proof`](super::Proof).
#[derive(
    CanonicalDeserialize,
    CanonicalSerialize,
    Clone,
    Debug,
    Default,
    Eq,
    PartialEq,
)]
pub struct ProofEvaluations<F: Field> {
    // Evaluation of the witness polynomial for the left wire at `z`
    pub(crate) a_eval: F,
    // Evaluation of the witness polynomial for the right wire at `z`
    pub(crate) b_eval: F,
    // Evaluation of the witness polynomial for the output wire at `z`
    pub(crate) c_eval: F,
    // Evaluation of the witness polynomial for the fourth wire at `z`
    pub(crate) d_eval: F,
    //
    pub(crate) a_next_eval: F,
    //
    pub(crate) b_next_eval: F,
    // Evaluation of the witness polynomial for the fourth wire at `z * root of
    // unity`
    pub(crate) d_next_eval: F,
    // Evaluation of the arithmetic selector polynomial at `z`
    pub(crate) q_arith_eval: F,
    //
    pub(crate) q_c_eval: F,
    //
    pub(crate) q_l_eval: F,
    //
    pub(crate) q_r_eval: F,
    // Evaluation of the left sigma polynomial at `z`
    pub(crate) left_sigma_eval: F,
    // Evaluation of the right sigma polynomial at `z`
    pub(crate) right_sigma_eval: F,
    // Evaluation of the out sigma polynomial at `z`
    pub(crate) out_sigma_eval: F,

    // Evaluation of the linearisation sigma polynomial at `z`
    pub(crate) lin_poly_eval: F,

    // (Shifted) Evaluation of the permutation polynomial at `z * root of
    // unity`
    pub(crate) perm_eval: F,
}

/// Compute the linearisation polynomial.
pub(crate) fn compute<
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
>(
    domain: &GeneralEvaluationDomain<E::Fr>,
    prover_key: &ProverKey<E::Fr, P>,
    alpha: &E::Fr,
    beta: &E::Fr,
    gamma: &E::Fr,
    range_separation_challenge: &E::Fr,
    logic_separation_challenge: &E::Fr,
    fixed_base_separation_challenge: &E::Fr,
    var_base_separation_challenge: &E::Fr,
    z_challenge: &E::Fr,
    w_l_poly: &DensePolynomial<E::Fr>,
    w_r_poly: &DensePolynomial<E::Fr>,
    w_o_poly: &DensePolynomial<E::Fr>,
    w_4_poly: &DensePolynomial<E::Fr>,
    t_x_poly: &DensePolynomial<E::Fr>,
    z_poly: &DensePolynomial<E::Fr>,
) -> (DensePolynomial<E::Fr>, Evaluations<E::Fr>) {
    // Compute evaluations
    let quot_eval = t_x_poly.evaluate(z_challenge);
    let a_eval = w_l_poly.evaluate(z_challenge);
    let b_eval = w_r_poly.evaluate(z_challenge);
    let c_eval = w_o_poly.evaluate(z_challenge);
    let d_eval = w_4_poly.evaluate(z_challenge);
    let left_sigma_eval =
        prover_key.permutation.left_sigma.0.evaluate(z_challenge);
    let right_sigma_eval =
        prover_key.permutation.right_sigma.0.evaluate(z_challenge);
    let out_sigma_eval =
        prover_key.permutation.out_sigma.0.evaluate(z_challenge);
    let q_arith_eval = prover_key.arithmetic.q_arith.0.evaluate(z_challenge);
    let q_c_eval = prover_key.constant_selector.0.evaluate(z_challenge);
    let q_l_eval = prover_key.left_selector.0.evaluate(z_challenge);
    let q_r_eval = prover_key.right_selector.0.evaluate(z_challenge);

    let group_gen = domain.group_gen();
    let a_next_eval = w_l_poly.evaluate(&(*z_challenge * group_gen));
    let b_next_eval = w_r_poly.evaluate(&(*z_challenge * group_gen));
    let d_next_eval = w_4_poly.evaluate(&(*z_challenge * group_gen));
    let perm_eval = z_poly.evaluate(&(*z_challenge * group_gen));

    let f_1 = compute_circuit_satisfiability::<E, P>(
        (
            range_separation_challenge,
            logic_separation_challenge,
            fixed_base_separation_challenge,
            var_base_separation_challenge,
        ),
        a_eval,
        b_eval,
        c_eval,
        d_eval,
        a_next_eval,
        b_next_eval,
        d_next_eval,
        q_arith_eval,
        q_c_eval,
        q_l_eval,
        q_r_eval,
        prover_key,
    );

    let f_2 = prover_key.permutation.compute_linearisation(
        *z_challenge,
        (*alpha, *beta, *gamma),
        (a_eval, b_eval, c_eval, d_eval),
        (left_sigma_eval, right_sigma_eval, out_sigma_eval),
        perm_eval,
        z_poly,
    );

    let lin_poly = &f_1 + &f_2;

    // Evaluate linearisation polynomial at z_challenge
    let lin_poly_eval = lin_poly.evaluate(z_challenge);

    (
        lin_poly,
        Evaluations {
            proof: ProofEvaluations {
                a_eval,
                b_eval,
                c_eval,
                d_eval,
                a_next_eval,
                b_next_eval,
                d_next_eval,
                q_arith_eval,
                q_c_eval,
                q_l_eval,
                q_r_eval,
                left_sigma_eval,
                right_sigma_eval,
                out_sigma_eval,
                lin_poly_eval,
                perm_eval,
            },
            quot_eval,
        },
    )
}

fn compute_circuit_satisfiability<
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
>(
    (
        range_separation_challenge,
        logic_separation_challenge,
        fixed_base_separation_challenge,
        var_base_separation_challenge,
    ): (&E::Fr, &E::Fr, &E::Fr, &E::Fr),
    a_eval: E::Fr,
    b_eval: E::Fr,
    c_eval: E::Fr,
    d_eval: E::Fr,
    a_next_eval: E::Fr,
    b_next_eval: E::Fr,
    d_next_eval: E::Fr,
    q_arith_eval: E::Fr,
    q_c_eval: E::Fr,
    q_l_eval: E::Fr,
    q_r_eval: E::Fr,
    prover_key: &ProverKey<E::Fr, P>,
) -> DensePolynomial<E::Fr> {
    let values = GateValues {
        left: a_eval,
        right: b_eval,
        output: c_eval,
        fourth: d_eval,
        left_next: a_next_eval,
        right_next: b_next_eval,
        fourth_next: d_next_eval,
        left_selector: q_l_eval,
        right_selector: q_r_eval,
        constant_selector: q_c_eval,
    };

    let a = prover_key.arithmetic.compute_linearisation(
        a_eval,
        b_eval,
        c_eval,
        d_eval,
        q_arith_eval,
    );

    let b = Range::linearisation_term(
        &prover_key.range_selector.0,
        *range_separation_challenge,
        values,
    );

    let c = Logic::linearisation_term(
        &prover_key.logic_selector.0,
        *logic_separation_challenge,
        values,
    );

    let d = FixedBaseScalarMul::<_, P>::linearisation_term(
        &prover_key.fixed_group_add_selector.0,
        *fixed_base_separation_challenge,
        values,
    );

    let e = CurveAddition::<_, P>::linearisation_term(
        &prover_key.variable_group_add_selector.0,
        *var_base_separation_challenge,
        values,
    );

    let mut linearisation_poly = &a + &b;
    linearisation_poly += &c;
    linearisation_poly += &d;
    linearisation_poly += &e;

    linearisation_poly
}
