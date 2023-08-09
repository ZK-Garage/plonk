// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use crate::{
    error::Error,
    proof_system::{
        ecc::{CurveAddition, FixedBaseScalarMul},
        logic::Logic,
        range::Range,
        widget::GateConstraint,
        ProverKey,
    },
};
use ark_ec::TEModelParameters;
use ark_ff::{FftField, PrimeField};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
    UVPolynomial,
};

use super::{
    ecc::{CAVals, FBSMVals},
    linearisation_poly::CustomEvaluations,
    logic::LogicVals,
    range::RangeVals,
    CustomValues, WitnessValues,
};

/// Computes the Quotient [`DensePolynomial`] given the [`EvaluationDomain`], a
/// [`ProverKey`], and some other info.
pub fn compute<F, P>(
    domain: &GeneralEvaluationDomain<F>,
    prover_key: &ProverKey<F>,
    z_poly: &DensePolynomial<F>,
    z2_poly: &DensePolynomial<F>,
    w_l_poly: &DensePolynomial<F>,
    w_r_poly: &DensePolynomial<F>,
    w_o_poly: &DensePolynomial<F>,
    w_4_poly: &DensePolynomial<F>,
    public_inputs_poly: &DensePolynomial<F>,
    f_poly: &DensePolynomial<F>,
    table_poly: &DensePolynomial<F>,
    h1_poly: &DensePolynomial<F>,
    h2_poly: &DensePolynomial<F>,
    alpha: &F,
    beta: &F,
    gamma: &F,
    delta: &F,
    epsilon: &F,
    zeta: &F,
    range_challenge: &F,
    logic_challenge: &F,
    fixed_base_challenge: &F,
    var_base_challenge: &F,
    lookup_challenge: &F,
) -> Result<DensePolynomial<F>, Error>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    let domain_8n = GeneralEvaluationDomain::<F>::new(8 * domain.size())
        .ok_or(Error::InvalidEvalDomainSize {
        log_size_of_group: (8 * domain.size()).trailing_zeros(),
        adicity:
            <<F as FftField>::FftParams as ark_ff::FftParameters>::TWO_ADICITY,
    })?;

    let l1_poly = compute_first_lagrange_poly_scaled(domain, F::one());
    let l1_eval_8n = domain_8n.coset_fft(&l1_poly);

    let mut z_eval_8n = domain_8n.coset_fft(z_poly);
    z_eval_8n.push(z_eval_8n[0]);
    z_eval_8n.push(z_eval_8n[1]);
    z_eval_8n.push(z_eval_8n[2]);
    z_eval_8n.push(z_eval_8n[3]);
    z_eval_8n.push(z_eval_8n[4]);
    z_eval_8n.push(z_eval_8n[5]);
    z_eval_8n.push(z_eval_8n[6]);
    z_eval_8n.push(z_eval_8n[7]);

    let mut wl_eval_8n = domain_8n.coset_fft(w_l_poly);
    wl_eval_8n.push(wl_eval_8n[0]);
    wl_eval_8n.push(wl_eval_8n[1]);
    wl_eval_8n.push(wl_eval_8n[2]);
    wl_eval_8n.push(wl_eval_8n[3]);
    wl_eval_8n.push(wl_eval_8n[4]);
    wl_eval_8n.push(wl_eval_8n[5]);
    wl_eval_8n.push(wl_eval_8n[6]);
    wl_eval_8n.push(wl_eval_8n[7]);

    let mut wr_eval_8n = domain_8n.coset_fft(w_r_poly);
    wr_eval_8n.push(wr_eval_8n[0]);
    wr_eval_8n.push(wr_eval_8n[1]);
    wr_eval_8n.push(wr_eval_8n[2]);
    wr_eval_8n.push(wr_eval_8n[3]);
    wr_eval_8n.push(wr_eval_8n[4]);
    wr_eval_8n.push(wr_eval_8n[5]);
    wr_eval_8n.push(wr_eval_8n[6]);
    wr_eval_8n.push(wr_eval_8n[7]);

    let wo_eval_8n = domain_8n.coset_fft(w_o_poly);

    let mut w4_eval_8n = domain_8n.coset_fft(w_4_poly);
    w4_eval_8n.push(w4_eval_8n[0]);
    w4_eval_8n.push(w4_eval_8n[1]);
    w4_eval_8n.push(w4_eval_8n[2]);
    w4_eval_8n.push(w4_eval_8n[3]);
    w4_eval_8n.push(w4_eval_8n[4]);
    w4_eval_8n.push(w4_eval_8n[5]);
    w4_eval_8n.push(w4_eval_8n[6]);
    w4_eval_8n.push(w4_eval_8n[7]);

    let mut z2_eval_8n = domain_8n.coset_fft(z2_poly);
    z2_eval_8n.push(z2_eval_8n[0]);
    z2_eval_8n.push(z2_eval_8n[1]);
    z2_eval_8n.push(z2_eval_8n[2]);
    z2_eval_8n.push(z2_eval_8n[3]);
    z2_eval_8n.push(z2_eval_8n[4]);
    z2_eval_8n.push(z2_eval_8n[5]);
    z2_eval_8n.push(z2_eval_8n[6]);
    z2_eval_8n.push(z2_eval_8n[7]);

    let f_eval_8n = domain_8n.coset_fft(f_poly);

    let mut table_eval_8n = domain_8n.coset_fft(table_poly);
    table_eval_8n.push(table_eval_8n[0]);
    table_eval_8n.push(table_eval_8n[1]);
    table_eval_8n.push(table_eval_8n[2]);
    table_eval_8n.push(table_eval_8n[3]);
    table_eval_8n.push(table_eval_8n[4]);
    table_eval_8n.push(table_eval_8n[5]);
    table_eval_8n.push(table_eval_8n[6]);
    table_eval_8n.push(table_eval_8n[7]);

    let mut h1_eval_8n = domain_8n.coset_fft(h1_poly);
    h1_eval_8n.push(h1_eval_8n[0]);
    h1_eval_8n.push(h1_eval_8n[1]);
    h1_eval_8n.push(h1_eval_8n[2]);
    h1_eval_8n.push(h1_eval_8n[3]);
    h1_eval_8n.push(h1_eval_8n[4]);
    h1_eval_8n.push(h1_eval_8n[5]);
    h1_eval_8n.push(h1_eval_8n[6]);
    h1_eval_8n.push(h1_eval_8n[7]);

    let h2_eval_8n = domain_8n.coset_fft(h2_poly);

    let gate_constraints = compute_gate_constraint_satisfiability::<F, P>(
        domain,
        *range_challenge,
        *logic_challenge,
        *fixed_base_challenge,
        *var_base_challenge,
        prover_key,
        &wl_eval_8n,
        &wr_eval_8n,
        &wo_eval_8n,
        &w4_eval_8n,
        public_inputs_poly,
    )?;

    let permutation = compute_permutation_checks::<F>(
        domain,
        prover_key,
        &wl_eval_8n,
        &wr_eval_8n,
        &wo_eval_8n,
        &w4_eval_8n,
        &z_eval_8n,
        *alpha,
        *beta,
        *gamma,
    )?;

    let lookup = prover_key.lookup.compute_lookup_quotient_term(
        domain,
        &wl_eval_8n,
        &wr_eval_8n,
        &wo_eval_8n,
        &w4_eval_8n,
        &f_eval_8n,
        &table_eval_8n,
        &h1_eval_8n,
        &h2_eval_8n,
        &z2_eval_8n,
        &l1_eval_8n,
        *delta,
        *epsilon,
        *zeta,
        *lookup_challenge,
    )?;

    let quotient = (0..domain_8n.size())
        .map(|i| {
            let numerator = gate_constraints[i] + permutation[i] + lookup[i];
            let denominator = prover_key.v_h_coset_8n()[i];
            numerator * denominator.inverse().unwrap()
        })
        .collect::<Vec<_>>();

    Ok(DensePolynomial::from_coefficients_vec(
        domain_8n.coset_ifft(&quotient),
    ))
}

/// Computes contribution to the quotient polynomial that ensures
/// the gate constraints are satisfied.
fn compute_gate_constraint_satisfiability<F, P>(
    domain: &GeneralEvaluationDomain<F>,
    range_challenge: F,
    logic_challenge: F,
    fixed_base_challenge: F,
    var_base_challenge: F,
    prover_key: &ProverKey<F>,
    wl_eval_8n: &[F],
    wr_eval_8n: &[F],
    wo_eval_8n: &[F],
    w4_eval_8n: &[F],
    pi_poly: &DensePolynomial<F>,
) -> Result<Vec<F>, Error>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    let domain_8n = GeneralEvaluationDomain::<F>::new(8 * domain.size())
        .ok_or(Error::InvalidEvalDomainSize {
        log_size_of_group: (8 * domain.size()).trailing_zeros(),
        adicity:
            <<F as FftField>::FftParams as ark_ff::FftParameters>::TWO_ADICITY,
    })?;
    let pi_eval_8n = domain_8n.coset_fft(pi_poly);

    // TODO Eliminate contribution of unused gates
    Ok((0..domain_8n.size())
        .map(|i| {
            let wit_vals = WitnessValues {
                a_val: wl_eval_8n[i],
                b_val: wr_eval_8n[i],
                c_val: wo_eval_8n[i],
                d_val: w4_eval_8n[i],
            };

            let custom_vals = CustomEvaluations {
                vals: vec![
                    ("a_next_eval".to_string(), wl_eval_8n[i + 8]),
                    ("b_next_eval".to_string(), wr_eval_8n[i + 8]),
                    ("d_next_eval".to_string(), w4_eval_8n[i + 8]),
                    ("q_l_eval".to_string(), prover_key.arithmetic.q_l.1[i]),
                    ("q_r_eval".to_string(), prover_key.arithmetic.q_r.1[i]),
                    ("q_c_eval".to_string(), prover_key.arithmetic.q_c.1[i]),
                    // prob don't need the follow but nonetheless...
                    ("q_hl_eval".to_string(), prover_key.arithmetic.q_hl.1[i]),
                    ("q_hr_eval".to_string(), prover_key.arithmetic.q_hr.1[i]),
                    ("q_h4_eval".to_string(), prover_key.arithmetic.q_hr.1[i]),
                ],
            };

            let arithmetic =
                prover_key.arithmetic.compute_quotient_i(i, wit_vals);

            let range = Range::quotient_term(
                prover_key.range_selector.1[i],
                range_challenge,
                wit_vals,
                RangeVals::from_evaluations(&custom_vals),
            );

            let logic = Logic::quotient_term(
                prover_key.logic_selector.1[i],
                logic_challenge,
                wit_vals,
                LogicVals::from_evaluations(&custom_vals),
            );

            let fixed_base_scalar_mul =
                FixedBaseScalarMul::<_, P>::quotient_term(
                    prover_key.fixed_group_add_selector.1[i],
                    fixed_base_challenge,
                    wit_vals,
                    FBSMVals::from_evaluations(&custom_vals),
                );

            let curve_addition = CurveAddition::<_, P>::quotient_term(
                prover_key.variable_group_add_selector.1[i],
                var_base_challenge,
                wit_vals,
                CAVals::from_evaluations(&custom_vals),
            );

            (arithmetic + pi_eval_8n[i])
                + range
                + logic
                + fixed_base_scalar_mul
                + curve_addition
        })
        .collect())
}

/// Computes the permutation contribution to the quotient polynomial over
/// `domain`.
fn compute_permutation_checks<F>(
    domain: &GeneralEvaluationDomain<F>,
    prover_key: &ProverKey<F>,
    wl_eval_8n: &[F],
    wr_eval_8n: &[F],
    wo_eval_8n: &[F],
    w4_eval_8n: &[F],
    z_eval_8n: &[F],
    alpha: F,
    beta: F,
    gamma: F,
) -> Result<Vec<F>, Error>
where
    F: PrimeField,
{
    let domain_8n = GeneralEvaluationDomain::<F>::new(8 * domain.size())
        .ok_or(Error::InvalidEvalDomainSize {
        log_size_of_group: (8 * domain.size()).trailing_zeros(),
        adicity:
            <<F as FftField>::FftParams as ark_ff::FftParameters>::TWO_ADICITY,
    })?;
    let l1_poly_alpha =
        compute_first_lagrange_poly_scaled(domain, alpha.square());
    let l1_alpha_sq_evals = domain_8n.coset_fft(&l1_poly_alpha.coeffs);
    Ok((0..domain_8n.size())
        .map(|i| {
            prover_key.permutation.compute_quotient_i(
                i,
                wl_eval_8n[i],
                wr_eval_8n[i],
                wo_eval_8n[i],
                w4_eval_8n[i],
                z_eval_8n[i],
                z_eval_8n[i + 8],
                alpha,
                l1_alpha_sq_evals[i],
                beta,
                gamma,
            )
        })
        .collect())
}

/// Computes the first lagrange polynomial with the given `scale` over `domain`.
fn compute_first_lagrange_poly_scaled<F>(
    domain: &GeneralEvaluationDomain<F>,
    scale: F,
) -> DensePolynomial<F>
where
    F: FftField,
{
    let mut x_evals = vec![F::zero(); domain.size()];
    x_evals[0] = scale;
    domain.ifft_in_place(&mut x_evals);
    DensePolynomial::from_coefficients_vec(x_evals)
}
