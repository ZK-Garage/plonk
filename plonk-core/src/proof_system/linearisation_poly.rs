// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::proof_system::ecc::{CurveAddition, FixedBaseScalarMul};
use crate::proof_system::logic::Logic;
use crate::proof_system::range::Range;
use crate::proof_system::widget::GateConstraint;
use crate::proof_system::GateValues;
use crate::util::EvaluationDomainExt;
use crate::proof_system::ProverKey;
use ark_ec::{PairingEngine, TEModelParameters};
use ark_ff::Field;
use ark_ff::PrimeField;
use ark_poly::EvaluationDomain;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
    Polynomial,
};
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write,
};
use crate::error::Error;
use ark_ec::TEModelParameters;
use ark_ff::{Field, PrimeField};

/// Polynomial Evaluations
///
/// This `struct` keeps track of polynomial evaluations at points `z` and/or `z
/// * w` where `w` is a root of unit.
pub struct Evaluations<F>
where
    F: Field,
{
    /// Proof-relevant Evaluations
    pub proof: ProofEvaluations<F>,

    /// Evaluation of the linearisation sigma polynomial at `z`.
    pub quot_eval: F,
}

/// Subset of all of the evaluations. These evaluations
/// are added to the [`Proof`](super::Proof).
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Default, Eq, PartialEq)]
pub struct ProofEvaluations<F>
where
    F: Field,
{
    /// Evaluation of the witness polynomial for the left wire at `z`.
    pub a_eval: F,

    /// Evaluation of the witness polynomial for the right wire at `z`.
    pub b_eval: F,

    /// Evaluation of the witness polynomial for the output wire at `z`.
    pub c_eval: F,

    /// Evaluation of the witness polynomial for the fourth wire at `z`.
    pub d_eval: F,

    /// Evaluation of the witness polynomial for the left wire at `z * omega`
    /// where `omega` is a root of unity.
    pub a_next_eval: F,

    /// Evaluation of the witness polynomial for the right wire at `z * omega`
    /// where `omega` is a root of unity.
    pub b_next_eval: F,

    /// Evaluation of the witness polynomial for the fourth wire at `z * omega`
    /// where `omega` is a root of unity.
    pub d_next_eval: F,

    /// Evaluation of the arithmetic selector polynomial at `z`.
    pub q_arith_eval: F,

    /// Evaluation of the constant selector polynomial at `z`.
    pub q_c_eval: F,

    /// Evaluation of the left selector polynomial at `z`.
    pub q_l_eval: F,

    /// Evaluation of the right selector polynomial at `z`.
    pub q_r_eval: F,

    /// Evaluation of the lookup selector polynomial at `z`.
    pub q_lookup_eval: F,

    /// Evaluation of the left sigma polynomial at `z`.
    pub left_sigma_eval: F,

    /// Evaluation of the right sigma polynomial at `z`.
    pub right_sigma_eval: F,

    /// Evaluation of the out sigma polynomial at `z`.
    pub out_sigma_eval: F,

    /// Evaluation of the linearisation sigma polynomial at `z`.
    pub linearisation_polynomial_eval: F,

    /// Evaluation of the permutation polynomial at `z * omega` where `omega`
    /// is a root of unity.
    pub permutation_eval: F,

    // (Shifted) Evaluation of the lookup permutation polynomial at `z * root
    // of unity`
    pub lookup_perm_eval: F,

    /// Evaluations of the first half of sorted plonkup poly at `z`
    pub h_1_eval: F,

    /// (Shifted) Evaluations of the first half of sorted plonkup poly at `z *
    /// root of unity`
    pub h_1_next_eval: F,

    /// (Shifted) Evaluations of the second half of sorted plonkup poly at `z *
    /// root of unity`
    pub h_2_eval: F,

    /// Evaluations of the query polynomial at `z`
    pub f_eval: F,

    /// Evaluations of the table polynomial at `z`
    pub table_eval: F,

    /// Evaluations of the table polynomial at `z * root of unity`
    pub table_next_eval: F,
}

/// Compute the linearisation polynomial.
pub fn compute<F, P>(
    domain: &GeneralEvaluationDomain<F>,
    prover_key: &ProverKey<F>,
    alpha: &F,
    beta: &F,
    gamma: &F,
    zeta: &F,
    range_separation_challenge: &F,
    logic_separation_challenge: &F,
    fixed_base_separation_challenge: &F,
    var_base_separation_challenge: &F,
    lookup_separation_challenge: &F,
    z_challenge: &F,
    w_l_poly: &DensePolynomial<F>,
    w_r_poly: &DensePolynomial<F>,
    w_o_poly: &DensePolynomial<F>,
    w_4_poly: &DensePolynomial<F>,
    quot_poly: &DensePolynomial<F>,
    z_poly: &DensePolynomial<F>,
    z_2_poly: &DensePolynomial<F>,
    f_poly: &DensePolynomial<F>,
    h_1_poly: &DensePolynomial<F>,
    h_2_poly: &DensePolynomial<F>,
    table_poly: &DensePolynomial<F>,
) -> Result<(DensePolynomial<F>, Evaluations<F>), Error>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    let quot_eval = quot_poly.evaluate(z_challenge);
    let a_eval = w_l_poly.evaluate(z_challenge);
    let b_eval = w_r_poly.evaluate(z_challenge);
    let c_eval = w_o_poly.evaluate(z_challenge);
    let d_eval = w_4_poly.evaluate(z_challenge);
    let f_eval = f_poly.evaluate(z_challenge);
    let left_sigma_eval =
        prover_key.permutation.left_sigma.0.evaluate(z_challenge);
    let right_sigma_eval =
        prover_key.permutation.right_sigma.0.evaluate(z_challenge);
    let out_sigma_eval =
        prover_key.permutation.out_sigma.0.evaluate(z_challenge);
    let q_arith_eval = prover_key.arithmetic.q_arith.0.evaluate(z_challenge);
    let q_c_eval = prover_key.arithmetic.q_c.0.evaluate(z_challenge);
    let q_l_eval = prover_key.arithmetic.q_l.0.evaluate(z_challenge);
    let q_r_eval = prover_key.arithmetic.q_r.0.evaluate(z_challenge);
    let q_lookup_eval = prover_key.lookup.q_lookup.0.evaluate(z_challenge);
    let h_1_eval = h_1_poly.evaluate(z_challenge);
    let h_2_eval = h_2_poly.evaluate(z_challenge);
    let table_eval = table_poly.evaluate(z_challenge);

    let omega = domain.group_gen();
    let shifted_z_challenge = *z_challenge * omega;

    let a_next_eval = w_l_poly.evaluate(&shifted_z_challenge);
    let b_next_eval = w_r_poly.evaluate(&shifted_z_challenge);
    let d_next_eval = w_4_poly.evaluate(&shifted_z_challenge);
    let permutation_eval = z_poly.evaluate(&shifted_z_challenge);
    let lookup_perm_eval = z_2_poly.evaluate(&shifted_z_challenge);
    let h_1_next_eval = h_1_poly.evaluate(&shifted_z_challenge);
    let table_next_eval = table_poly.evaluate(&shifted_z_challenge);

    let gate_constraints = compute_gate_constraint_satisfiability::<F, P>(
        range_separation_challenge,
        logic_separation_challenge,
        fixed_base_separation_challenge,
        var_base_separation_challenge,
        logic_separation_challenge,
        a_eval,
        b_eval,
        c_eval,
        d_eval,
        f_eval,
        a_next_eval,
        b_next_eval,
        d_next_eval,
        q_arith_eval,
        q_c_eval,
        q_l_eval,
        q_r_eval,
        zeta,
        prover_key,
    );

    let permutation = prover_key.permutation.compute_linearisation(
        domain.size(),
        *z_challenge,
        (*alpha, *beta, *gamma),
        (a_eval, b_eval, c_eval, d_eval),
        (left_sigma_eval, right_sigma_eval, out_sigma_eval),
        permutation_eval,
        z_poly,
    )?;

    let linearisation_polynomial = gate_constraints + permutation;
    let linearisation_polynomial_eval =
        linearisation_polynomial.evaluate(z_challenge);

    Ok((
        linearisation_polynomial,
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
                q_lookup_eval,
                left_sigma_eval,
                right_sigma_eval,
                out_sigma_eval,
                linearisation_polynomial_eval,
                permutation_eval,
                lookup_perm_eval,
                h_1_eval,
                h_1_next_eval,
                h_2_eval,
                f_eval,
                table_eval,
                table_next_eval,
                
            },
            quot_eval,
        },
    ))
}

/// Computes the gate constraint satisfiability portion of the linearisation
/// polynomial.
fn compute_gate_constraint_satisfiability<F, P>(
    range_separation_challenge: &F,
    logic_separation_challenge: &F,
    fixed_base_separation_challenge: &F,
    var_base_separation_challenge: &F,
    lookup_separation_challenge: &F,
    a_eval: F,
    b_eval: F,
    c_eval: F,
    d_eval: F,
    f_eval: F,
    a_next_eval: F,
    b_next_eval: F,
    d_next_eval: F,
    q_arith_eval: F,
    q_c_eval: F,
    q_l_eval: F,
    q_r_eval: F,
    zeta: &F,
    prover_key: &ProverKey<F>,
) -> DensePolynomial<F>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
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

    let arithmetic = prover_key.arithmetic.compute_linearisation(
        a_eval,
        b_eval,
        c_eval,
        d_eval,
        q_arith_eval,
    );
    let range = Range::linearisation_term(
        &prover_key.range_selector.0,
        *range_separation_challenge,
        values,
    );

    let logic = Logic::linearisation_term(
        &prover_key.logic_selector.0,
        *logic_separation_challenge,
        values,
    );

    let fixed_base_scalar_mul = FixedBaseScalarMul::<F, P>::linearisation_term(
        &prover_key.fixed_group_add_selector.0,
        *fixed_base_separation_challenge,
        values,
    );

    let curve_addition = CurveAddition::<F, P>::linearisation_term(
        &prover_key.variable_group_add_selector.0,
        *var_base_separation_challenge,
        values,
    );

    let lookup = prover_key.lookup.compute_linearization(a_eval, b_eval, c_eval, d_eval, f_eval, *zeta, *lookup_separation_challenge);

    arithmetic + range + logic + fixed_base_scalar_mul + curve_addition + lookup
}
