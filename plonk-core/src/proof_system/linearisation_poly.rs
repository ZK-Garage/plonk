// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::{
    error::Error,
    label_eval,
    parameters::CircuitParameters,
    proof_system::{
        logic::{Logic, LogicVals},
        proof,
        range::{Range, RangeVals},
        widget::GateConstraint,
        CustomValues, ProverKey, WitnessValues,
    },
    util::EvaluationDomainExt,
};
use ark_ff::{Field, One};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
    Polynomial,
};
use ark_serialize::{
    CanonicalDeserialize, CanonicalSerialize, Read, SerializationError, Write,
};

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

// Probably all of these should go into CustomEvals
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Default, Eq, PartialEq)]
pub struct LookupEvaluations<F>
where
    F: Field,
{
    pub q_lookup_eval: F,
    // (Shifted) Evaluation of the lookup permutation polynomial at `z * root
    // of unity`
    pub z2_next_eval: F,

    /// Evaluations of the first half of sorted plonkup poly at `z`
    pub h1_eval: F,

    /// (Shifted) Evaluations of the even indexed half of sorted plonkup poly
    /// at `z root of unity
    pub h1_next_eval: F,

    /// Evaluations of the odd indexed half of sorted plonkup poly at `z
    /// root of unity
    pub h2_eval: F,

    /// Evaluations of the query polynomial at `z`
    pub f_eval: F,

    /// Evaluations of the table polynomial at `z`
    pub table_eval: F,

    /// Evaluations of the table polynomial at `z * root of unity`
    pub table_next_eval: F,
}

/// Subset of the [`ProofEvaluations`]. Evaluations at `z`  or `z *w` where `w`
/// is the nth root of unity of selectors polynomials needed for custom gates
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Default, Eq, PartialEq)]
pub struct CustomEvaluations<F>
where
    F: Field,
{
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
    pub fn add(&mut self, label: &str, poly: DensePolynomial<F>, point: F) {
        if let Some(_l) = &self.vals.iter().find(|entry| entry.0 == label) {
        } else {
            let eval = poly.evaluate(&point);
            let _ = &self.vals.push((label.to_string(), eval));
        }
    }
}

/// Set of evaluations that form the [`Proof`](super::Proof).
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

    /// Lookup evaluations
    pub lookup_evals: LookupEvaluations<F>,

    /// Evaluations needed for custom gates. This includes selector polynomials
    /// and evaluations of wire polynomials at an offset
    pub custom_evals: CustomEvaluations<F>,
}

/// Compute the linearisation polynomial.
pub fn compute<P>(
    domain: &GeneralEvaluationDomain<P::ScalarField>,
    prover_key: &ProverKey<P::ScalarField>,
    alpha: &P::ScalarField,
    beta: &P::ScalarField,
    gamma: &P::ScalarField,
    delta: &P::ScalarField,
    epsilon: &P::ScalarField,
    zeta: &P::ScalarField,
    range_separation_challenge: &P::ScalarField,
    logic_separation_challenge: &P::ScalarField,
    fixed_base_separation_challenge: &P::ScalarField,
    var_base_separation_challenge: &P::ScalarField,
    lookup_separation_challenge: &P::ScalarField,
    z_challenge: &P::ScalarField,
    w_l_poly: &DensePolynomial<P::ScalarField>,
    w_r_poly: &DensePolynomial<P::ScalarField>,
    w_o_poly: &DensePolynomial<P::ScalarField>,
    w_4_poly: &DensePolynomial<P::ScalarField>,
    t_1_poly: &DensePolynomial<P::ScalarField>,
    t_2_poly: &DensePolynomial<P::ScalarField>,
    t_3_poly: &DensePolynomial<P::ScalarField>,
    t_4_poly: &DensePolynomial<P::ScalarField>,
    z_poly: &DensePolynomial<P::ScalarField>,
    z2_poly: &DensePolynomial<P::ScalarField>,
    f_poly: &DensePolynomial<P::ScalarField>,
    h1_poly: &DensePolynomial<P::ScalarField>,
    h2_poly: &DensePolynomial<P::ScalarField>,
    table_poly: &DensePolynomial<P::ScalarField>,
) -> Result<
    (
        DensePolynomial<P::ScalarField>,
        ProofEvaluations<P::ScalarField>,
    ),
    Error,
>
where
    P: CircuitParameters,
{
    let n = domain.size();
    let omega = domain.group_gen();
    let shifted_z_challenge = *z_challenge * omega;

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

    // Lookup selector evaluation
    let q_lookup_eval = prover_key.lookup.q_lookup.0.evaluate(z_challenge);

    // Custom gate evaluations
    let q_c_eval = prover_key.arithmetic.q_c.0.evaluate(z_challenge);
    let q_l_eval = prover_key.arithmetic.q_l.0.evaluate(z_challenge);
    let q_r_eval = prover_key.arithmetic.q_r.0.evaluate(z_challenge);
    let a_next_eval = w_l_poly.evaluate(&shifted_z_challenge);
    let b_next_eval = w_r_poly.evaluate(&shifted_z_challenge);
    let d_next_eval = w_4_poly.evaluate(&shifted_z_challenge);

    let custom_evals = CustomEvaluations {
        vals: vec![
            label_eval!(q_arith_eval),
            label_eval!(q_c_eval),
            label_eval!(q_l_eval),
            label_eval!(q_r_eval),
            label_eval!(a_next_eval),
            label_eval!(b_next_eval),
            label_eval!(d_next_eval),
        ],
    };

    let z2_next_eval = z2_poly.evaluate(&shifted_z_challenge);
    let h1_eval = h1_poly.evaluate(z_challenge);
    let h1_next_eval = h1_poly.evaluate(&shifted_z_challenge);
    let h2_eval = h2_poly.evaluate(z_challenge);
    let f_eval = f_poly.evaluate(z_challenge);
    let table_eval = table_poly.evaluate(z_challenge);
    let table_next_eval = table_poly.evaluate(&shifted_z_challenge);

    // Compute the last term in the linearisation polynomial
    // (negative_quotient_term):
    // - Z_h(z_challenge) * [t_1(X) + z_challenge^n * t_2(X) + z_challenge^2n *
    //   t_3(X) + z_challenge^3n * t_4(X)]
    let vanishing_poly_eval =
        domain.evaluate_vanishing_polynomial(*z_challenge);
    let z_challenge_to_n = vanishing_poly_eval + P::ScalarField::one();
    let l1_eval = proof::compute_first_lagrange_evaluation(
        domain,
        &vanishing_poly_eval,
        z_challenge,
    );

    let lookup_evals = LookupEvaluations {
        q_lookup_eval,
        z2_next_eval,
        h1_eval,
        h1_next_eval,
        h2_eval,
        f_eval,
        table_eval,
        table_next_eval,
    };

    let gate_constraints = compute_gate_constraint_satisfiability::<P>(
        range_separation_challenge,
        logic_separation_challenge,
        fixed_base_separation_challenge,
        var_base_separation_challenge,
        &wire_evals,
        q_arith_eval,
        &custom_evals,
        prover_key,
    );

    let lookup = prover_key.lookup.compute_linearisation(
        l1_eval,
        a_eval,
        b_eval,
        c_eval,
        d_eval,
        f_eval,
        table_eval,
        table_next_eval,
        h1_next_eval,
        h2_eval,
        z2_next_eval,
        *delta,
        *epsilon,
        *zeta,
        z2_poly,
        h1_poly,
        *lookup_separation_challenge,
    );

    let permutation = prover_key.permutation.compute_linearisation(
        n,
        *z_challenge,
        (*alpha, *beta, *gamma),
        (a_eval, b_eval, c_eval, d_eval),
        (left_sigma_eval, right_sigma_eval, out_sigma_eval),
        permutation_eval,
        z_poly,
    )?;

    let quotient_term = &(&(&(&(&(&(t_4_poly * z_challenge_to_n)
        + t_3_poly)
        * z_challenge_to_n)
        + t_2_poly)
        * z_challenge_to_n)
        + t_1_poly)
        * vanishing_poly_eval;
    let negative_quotient_term = &quotient_term * (-P::ScalarField::one());

    let linearisation_polynomial =
        gate_constraints + permutation + lookup + negative_quotient_term;

    Ok((
        linearisation_polynomial,
        ProofEvaluations {
            wire_evals,
            perm_evals,
            lookup_evals,
            custom_evals,
        },
    ))
}

/// Computes the gate constraint satisfiability portion of the linearisation
/// polynomial.
fn compute_gate_constraint_satisfiability<P>(
    range_separation_challenge: &P::ScalarField,
    logic_separation_challenge: &P::ScalarField,
    fixed_base_separation_challenge: &P::ScalarField,
    var_base_separation_challenge: &P::ScalarField,
    wire_evals: &WireEvaluations<P::ScalarField>,
    q_arith_eval: P::ScalarField,
    custom_evals: &CustomEvaluations<P::ScalarField>,
    prover_key: &ProverKey<P::ScalarField>,
) -> DensePolynomial<P::ScalarField>
where
    P: CircuitParameters,
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
        RangeVals::from_evaluations(custom_evals),
    );

    let logic = Logic::linearisation_term(
        &prover_key.logic_selector.0,
        *logic_separation_challenge,
        wit_vals,
        LogicVals::from_evaluations(custom_evals),
    );

    let fbsm_custom_vals =
        <<P as CircuitParameters>::FixedBaseScalarMul as GateConstraint<
            P::ScalarField,
        >>::CustomVals::from_evaluations(custom_evals); 

    let fixed_base_scalar_mul = P::FixedBaseScalarMul::linearisation_term(
        &prover_key.fixed_group_add_selector.0,
        *fixed_base_separation_challenge,
        wit_vals,
        fbsm_custom_vals,
    );

    let ca_custom_vals =
        <<P as CircuitParameters>::CurveAddition as GateConstraint<
            P::ScalarField,
        >>::CustomVals::from_evaluations(custom_evals); 

    let curve_addition = P::CurveAddition::linearisation_term(
        &prover_key.variable_group_add_selector.0,
        *var_base_separation_challenge,
        wit_vals,
        ca_custom_vals,
    );

    arithmetic + range + logic + fixed_base_scalar_mul + curve_addition
}
