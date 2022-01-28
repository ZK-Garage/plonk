// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::{
    error::Error,
    label_eval,
    proof_system::{
        ecc::{CurveAddition, FixedBaseScalarMul},
        logic::Logic,
        range::Range,
        widget::GateConstraint,
        ProverKey, WitnessValues,
    },
    util::EvaluationDomainExt,
};
use ark_ec::TEModelParameters;
use ark_ff::{Field, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
    Polynomial,
};
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write,
};

use super::{
    ecc::{CAVals, FBSMVals},
    logic::LogicVals,
    range::RangeVals,
    CustomValues,
};

/// Polynomial Evaluations
///
/// This `struct` keeps track of polynomial evaluations at points `z` and/or `z
/// * w` where `w` is a root of unit.
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Default, Eq, PartialEq)]
pub struct Evaluations<F>
where
    F: Field,
{
    /// Proof-relevant Evaluations
    pub proof: ProofEvaluations<F>,

    /// Evaluation of the linearisation sigma polynomial at `z`.
    pub quot_eval: F,
}

/// Subset of the [`ProofEvaluations`]. Evaluations at `z` of the
/// wire polynomials
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct WireEvaluations<F>
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
}

/// Subset of the [`ProofEvaluations`]. Evaluations of the sigma and permutation
/// polynomials at `z`  or `z *w` where `w` is the nth root of unity.
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct PermutationEvaluations<F>
where
    F: Field,
{
    /// Evaluation of the left sigma polynomial at `z`.
    pub left_sigma_eval: F,

    /// Evaluation of the right sigma polynomial at `z`.
    pub right_sigma_eval: F,

    /// Evaluation of the out sigma polynomial at `z`.
    pub out_sigma_eval: F,

    /// Evaluation of the permutation polynomial at `z * omega` where `omega`
    /// is a root of unity.
    pub permutation_eval: F,
}

/// Subset of the [`ProofEvaluations`]. Evaluations at `z`  or `z *w` where `w`
/// is the nth root of unity of selectors polynomials needed for custom gates
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Default, Eq, PartialEq)]
pub struct CustomEvaluations<F>
where
    F: Field,
{
    // TODO: Decide best data structure
    // pub vals: HashMap<Into<String>, F>,
    pub vals: Vec<(String, F)>,
}

impl<F> CustomEvaluations<F>
where
    F: Field,
{
    /// Get the evaluation of the specified label.
    /// This funtions panics if the requested label is not found
    pub fn get(&self, label: &str) -> F {
        if let Some(result) = &self.vals.iter().find(|entry| entry.0 == label) {
            result.1
        } else {
            panic!("{} label not found in evaluations set", label)
        }
    }

    /// Add evaluation of poly at point if the label is not already
    /// in the set of evaluations
    pub fn add(label: String, poly: DensePolynomial<F>, point: F) {
        // TODO
        unimplemented!()
    }
}

/// Subset of all of the evaluations. These evaluations
/// are added to the [`Proof`](super::Proof).
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Default, Eq, PartialEq)]
pub struct ProofEvaluations<F>
where
    F: Field,
{
    /// Wire evaluations
    pub wire_evals: WireEvaluations<F>,

    /// Permutation and sigma polynomials evaluations
    pub perm_evals: PermutationEvaluations<F>,

    /// Evaluations needed for custom gates. This includes selector polynomials
    /// and evaluations of wire polynomials at an offset
    pub custom_evals: CustomEvaluations<F>,

    // The arithmetic gates is not considered custom
    // TODO It can be considered optional -> It should be removed when there
    // are no other gates (a purely arithmetic circuit)
    /// Evaluation of the arithmetic selector polynomial at `z`.
    pub q_arith_eval: F,

    // TODO This evaluation can be removed, making the adjustment for the
    // linearisation polynomial proposed in the last version of the print
    /// Evaluation of the linearisation sigma polynomial at `z`.
    pub linearisation_polynomial_eval: F,
}

/// Compute the linearisation polynomial.
pub fn compute<F, P>(
    domain: &GeneralEvaluationDomain<F>,
    prover_key: &ProverKey<F>,
    alpha: &F,
    beta: &F,
    gamma: &F,
    range_separation_challenge: &F,
    logic_separation_challenge: &F,
    fixed_base_separation_challenge: &F,
    var_base_separation_challenge: &F,
    z_challenge: &F,
    w_l_poly: &DensePolynomial<F>,
    w_r_poly: &DensePolynomial<F>,
    w_o_poly: &DensePolynomial<F>,
    w_4_poly: &DensePolynomial<F>,
    t_x_poly: &DensePolynomial<F>,
    z_poly: &DensePolynomial<F>,
) -> Result<(DensePolynomial<F>, Evaluations<F>), Error>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    let omega = domain.group_gen();
    let shifted_z_challenge = *z_challenge * omega;

    let quot_eval = t_x_poly.evaluate(z_challenge);

    // Wire evaluations
    let a_eval = w_l_poly.evaluate(z_challenge);
    let b_eval = w_r_poly.evaluate(z_challenge);
    let c_eval = w_o_poly.evaluate(z_challenge);
    let d_eval = w_4_poly.evaluate(z_challenge);
    let wire_evals = WireEvaluations {
        a_eval,
        b_eval,
        c_eval,
        d_eval,
    };

    // Permutation evaluations
    let left_sigma_eval =
        prover_key.permutation.left_sigma.0.evaluate(z_challenge);
    let right_sigma_eval =
        prover_key.permutation.right_sigma.0.evaluate(z_challenge);
    let out_sigma_eval =
        prover_key.permutation.out_sigma.0.evaluate(z_challenge);
    let permutation_eval = z_poly.evaluate(&shifted_z_challenge);
    let perm_evals = PermutationEvaluations {
        left_sigma_eval,
        right_sigma_eval,
        out_sigma_eval,
        permutation_eval,
    };

    // Arith selector evaluation
    let q_arith_eval = prover_key.arithmetic.q_arith.0.evaluate(z_challenge);

    // Custom gate evaluations
    let q_c_eval = prover_key.arithmetic.q_c.0.evaluate(z_challenge);
    let q_l_eval = prover_key.arithmetic.q_l.0.evaluate(z_challenge);
    let q_r_eval = prover_key.arithmetic.q_r.0.evaluate(z_challenge);

    let a_next_eval = w_l_poly.evaluate(&shifted_z_challenge);
    let b_next_eval = w_r_poly.evaluate(&shifted_z_challenge);
    let d_next_eval = w_4_poly.evaluate(&shifted_z_challenge);

    let custom_evals = CustomEvaluations {
        vals: vec![
            label_eval!(q_c_eval),
            label_eval!(q_l_eval),
            label_eval!(q_r_eval),
            label_eval!(a_next_eval),
            label_eval!(b_next_eval),
            label_eval!(d_next_eval),
        ],
    };

    let gate_constraints = compute_gate_constraint_satisfiability::<F, P>(
        range_separation_challenge,
        logic_separation_challenge,
        fixed_base_separation_challenge,
        var_base_separation_challenge,
        &wire_evals,
        q_arith_eval,
        &custom_evals,
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
                wire_evals,
                q_arith_eval,
                perm_evals,
                custom_evals,
                linearisation_polynomial_eval,
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
    wire_evals: &WireEvaluations<F>,
    q_arith_eval: F,
    custom_evals: &CustomEvaluations<F>,
    prover_key: &ProverKey<F>,
) -> DensePolynomial<F>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    let wit_vals = WitnessValues {
        a_val: wire_evals.a_eval,
        b_val: wire_evals.b_eval,
        c_val: wire_evals.c_eval,
        d_val: wire_evals.d_eval,
    };

    let arithmetic = prover_key.arithmetic.compute_linearisation(
        wire_evals.a_eval,
        wire_evals.b_eval,
        wire_evals.c_eval,
        wire_evals.d_eval,
        q_arith_eval,
    );

    let range = Range::linearisation_term(
        &prover_key.range_selector.0,
        *range_separation_challenge,
        wit_vals,
        RangeVals::from_evaluations(&custom_evals),
    );

    let logic = Logic::linearisation_term(
        &prover_key.logic_selector.0,
        *logic_separation_challenge,
        wit_vals,
        LogicVals::from_evaluations(&custom_evals),
    );

    let fixed_base_scalar_mul = FixedBaseScalarMul::<F, P>::linearisation_term(
        &prover_key.fixed_group_add_selector.0,
        *fixed_base_separation_challenge,
        wit_vals,
        FBSMVals::from_evaluations(&custom_evals),
    );

    let curve_addition = CurveAddition::<F, P>::linearisation_term(
        &prover_key.variable_group_add_selector.0,
        *var_base_separation_challenge,
        wit_vals,
        CAVals::from_evaluations(&custom_evals),
    );

    arithmetic + range + logic + fixed_base_scalar_mul + curve_addition
}
