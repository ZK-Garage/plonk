// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

use alloc::vec::Vec;
use ark_ec::{ModelParameters, TEModelParameters};
use ark_ff::{BigInteger, FftField, Field, FpParameters, PrimeField};
use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, univariate::DensePolynomial};
use hashbrown::HashMap;

use crate::{constraint_system::{StandardComposer, Variable}, lookup::LookupTable};

/// Returns an iterator over increasing powers of the given `scalar` starting
/// at `0`.
#[inline]
pub fn powers_of<F>(scalar: F) -> impl Iterator<Item = F>
where
    F: Field,
{
    core::iter::successors(Some(F::one()), move |p| Some(*p * scalar))
}

/// Evaluation Domain Extension Trait
pub trait EvaluationDomainExt<F>: EvaluationDomain<F>
where
    F: FftField,
{
    /// Returns the value of `log_2(self.size)`.
    fn log_size_of_group(&self) -> u32;

    /// Returns the inverse of the size in the field.
    fn size_inv(&self) -> F;

    /// Returns a fixed generator of the subgroup.
    fn group_gen(&self) -> F;

    /// Returns the inverse of the fixed generator of the subgroup.
    fn group_gen_inv(&self) -> F;

    /// Returns a fixed multiplicative generator of the finite field.
    fn generator_inv(&self) -> F;
}

impl<F> EvaluationDomainExt<F> for GeneralEvaluationDomain<F>
where
    F: FftField,
{
    #[inline]
    fn log_size_of_group(&self) -> u32 {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.log_size_of_group,
            GeneralEvaluationDomain::MixedRadix(domain) => {
                domain.log_size_of_group
            }
        }
    }

    #[inline]
    fn size_inv(&self) -> F {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.size_inv,
            GeneralEvaluationDomain::MixedRadix(domain) => domain.size_inv,
        }
    }

    #[inline]
    fn group_gen(&self) -> F {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.group_gen,
            GeneralEvaluationDomain::MixedRadix(domain) => domain.group_gen,
        }
    }

    #[inline]
    fn group_gen_inv(&self) -> F {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.group_gen_inv,
            GeneralEvaluationDomain::MixedRadix(domain) => domain.group_gen_inv,
        }
    }

    #[inline]
    fn generator_inv(&self) -> F {
        match self {
            GeneralEvaluationDomain::Radix2(domain) => domain.generator_inv,
            GeneralEvaluationDomain::MixedRadix(domain) => domain.generator_inv,
        }
    }
}

/// Get a pairing friendly curve scalar `E::Fr` from a scalar of the embedded
/// curve. Panics if the embedded scalar is greater than the modulus of the
/// pairing firendly curve scalar field
#[allow(dead_code)]
pub fn from_embedded_curve_scalar<F, P>(
    embedded_scalar: <P as ModelParameters>::ScalarField,
) -> F
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    let scalar_repr = embedded_scalar.into_repr();
    let modulus = <<F as PrimeField>::Params as FpParameters>::MODULUS;
    if modulus.num_bits() >= scalar_repr.num_bits() {
        let s = <<F as PrimeField>::BigInt as BigInteger>::from_bits_le(
            &scalar_repr.to_bits_le(),
        );
        assert!(s < modulus,
            "The embedded scalar exceeds the capacity representation of the outter curve scalar");
    } else {
        let m = <<P::ScalarField as PrimeField>::BigInt as BigInteger>::from_bits_le(
            &modulus.to_bits_le(),
        );
        assert!(scalar_repr < m,
            "The embedded scalar exceeds the capacity representation of the outter curve scalar");
    }
    F::from_le_bytes_mod_order(&scalar_repr.to_bytes_le())
}

/// Get a embedded curve scalar `P::ScalarField` from a scalar of the pariring
/// friendly curve. Panics if the pairing frindly curve scalar is greater than
/// the modulus of the embedded curve scalar field
#[allow(dead_code)]
pub(crate) fn to_embedded_curve_scalar<F, P>(pfc_scalar: F) -> P::ScalarField
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    let scalar_repr = pfc_scalar.into_repr();
    let modulus =
        <<P::ScalarField as PrimeField>::Params as FpParameters>::MODULUS;
    if modulus.num_bits() >= scalar_repr.num_bits() {
        let s = <<P::ScalarField as PrimeField>::BigInt as BigInteger>::from_bits_le(
            &scalar_repr.to_bits_le(),
        );
        assert!(s < modulus,
            "The embedded scalar exceeds the capacity representation of the outter curve scalar");
    } else {
        let m = <<F as PrimeField>::BigInt as BigInteger>::from_bits_le(
            &modulus.to_bits_le(),
        );
        assert!(scalar_repr < m,
            "The embedded scalar exceeds the capacity representation of the outter curve scalar");
    }
    P::ScalarField::from_le_bytes_mod_order(&scalar_repr.to_bytes_le())
}

/// Linear combination of a vector of values
///
/// For values [v_0, v_1,... v_k] returns:
/// v_0 + challenge * v_1 + ... + challenge^k  * v_k
pub fn lc<F>(values: Vec<F>, challenge: F) -> F
where
    F: Field,
{
    // Ensure valid challenge
    assert_ne!(challenge, F::zero());
    assert_ne!(challenge, F::one());

    values
        .iter()
        .rev()
        .fold(F::zero(), |acc, val| acc * challenge + *val)
}

/// Macro to quickly label polynomials
#[macro_export]
macro_rules! label_polynomial {
    ($poly:expr) => {
        ark_poly_commit::LabeledPolynomial::new(
            stringify!($poly).to_owned(),
            $poly.clone(),
            None,
            None,
        )
    };
}

/// Macro to quickly label polynomial commitments
#[macro_export]
macro_rules! label_commitment {
    ($comm:expr) => {
        ark_poly_commit::LabeledCommitment::new(
            stringify!($comm).to_owned(),
            $comm.clone(),
            None,
        )
    };
}

/// Macro to quickly label evaluations
#[macro_export]
macro_rules! label_eval {
    ($eval:expr) => {
        (stringify!($eval).to_owned(), $eval)
    };
}

/// Macro to get appropirate label
#[macro_export]
macro_rules! get_label {
    ($eval:expr) => {
        stringify!($comm).to_owned()
    };
}

/// abbreviates a field element to the first and last three hex digits
/// joined by ".."
pub fn abbreviate_field<F: Field>(f: &F) -> String {
    let whole_num_string = f
        .to_string()
        .split(|c| c == '(' || c == ')')
        .collect::<Vec<_>>()[1]
        .to_string();
    let len = whole_num_string.len();
    [&whole_num_string[..3], &whole_num_string[(len - 3)..]].join("##")
}

/// prints a vector of abbreviated field elements
pub fn abbreviate_vec<F: Field>(v: &Vec<F>) -> String {
    format!(
        "[{}]",
        v.iter()
            .map(|f| abbreviate_field(f))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn var_map_debug_print<F: Field>(hm: &HashMap<Variable, F>) -> String
where
    F: PrimeField,
{
    let mut var_indices =
        hm.keys().map(|&Variable(i)| i).collect::<Vec<_>>();
    var_indices.sort();

    var_indices
        .into_iter()
        .map(|i| {
            format!(
                "    {:<16} {}",
                format!("{:?}:", Variable(i)),
                abbreviate_field(
                    hm.get(&Variable(i)).unwrap()
                )
            )
        })
        .collect::<Vec<_>>()
        .join("\n")
}

/// prints polynomial form with abbreviated field elements
pub fn poly_coeff_debug_print<F: Field>(poly: &DensePolynomial<F>) -> String
where
    F: PrimeField,
{
    poly.coeffs
        .iter()
        .zip(0..poly.coeffs.len())
        .map(|(s, i)| {
            format!(
                "{} * x^{}",
                abbreviate_field(s),
                i
            )
        })
        .collect::<Vec<_>>()
        .join(" + ")
}

    /// prints polynomial evals with abbreviated field elements
    pub fn poly_eval_debug_print<F: Field>(poly: &DensePolynomial<F>, domain: GeneralEvaluationDomain<F>) -> String
    where
        F: PrimeField,
    {
        poly.evaluate_over_domain_by_ref(domain).evals
            .iter()
            .map(|f| abbreviate_field(f))
            .collect::<Vec<_>>()
            .join(", ")
    }

    /// prints both standard and eval forms
    pub fn poly_debug_print<F: PrimeField>(poly: &DensePolynomial<F>, domain: GeneralEvaluationDomain<F>)->String{
        format!("\t{}\n\t{}", poly_coeff_debug_print(poly), poly_eval_debug_print(poly, domain))
    }

/// prints wire variables in a table
pub fn wires_debug_print<F: Field>(wires: Vec<&Vec<Variable>>) -> String
where
    F: PrimeField,
{
    let wire_len = wires[0].len();

    let entires = (0..wire_len)
        .map(|i| {
            wires
                .iter()
                .map(|w| format!("{:<14}", format!("{:?}", w[i])))
                .collect::<Vec<_>>()
                .join("")
        })
        .collect::<Vec<_>>()
        .join("\n    ");

    let headings = format!(
        "{:<4}{:<14}{:<14}{:<14}{:<14}",
        "", "w_l", "w_r", "w_o", "w_4"
    );
    let debug_str = format!("{}\n    {}", headings, entires);

    debug_str
}

/// prints wire variables in a table
pub fn wires_values_debug_print<F: Field>(wires: Vec<&Vec<Variable>>, hm: &HashMap<Variable, F>) -> String
where
    F: PrimeField,
{
    let wire_len = wires[0].len();

    let entires = (0..wire_len)
        .map(|i| {
            wires
                .iter()
                .map(|w| format!("{:<14}", format!("{}", abbreviate_field(hm.get(&w[i]).unwrap()))))
                .collect::<Vec<_>>()
                .join("")
        })
        .collect::<Vec<_>>()
        .join("\n    ");

    let headings = format!(
        "{:<4}{:<14}{:<14}{:<14}{:<14}",
        "", "w_l", "w_r", "w_o", "w_4"
    );
    let debug_str = format!("{}\n    {}", headings, entires);

    debug_str
}

/// prints out composer in nicer form
pub fn debug_print<F: PrimeField, P: ModelParameters<BaseField = F>>(composer: StandardComposer<F, P>) -> String {
    let debug_str = format!(
        "\
        \n\
        circuit size: {}\n\
        {:<27} {}\n\
        {:<27} {}\n\
        {:<27} {}\n\
        {:<27} {}\n\
        {:<27} {}\n\
        {:<27} {}\n\
        {:<27} {}\n\
        {:<27} {}\n\
        {:<27} {}\n\
        {:<27} {}\n\
        {:<27} {}\n\
        {:<27} {}\n\n\
        public_inputs: {:?}\n\n\
        {}\n{}\n\n\
        {}\n{}\n\n\
        {}\n{}\n\n\
        {}\n{}\n\n\
        permutation: {:?}\n\n\
    ",
        composer.n,
        "q_m evals:",
        abbreviate_vec(&composer.q_m),
        "q_l evals:",
        abbreviate_vec(&composer.q_l),
        "q_r evals:",
        abbreviate_vec(&composer.q_r),
        "q_o evals:",
        abbreviate_vec(&composer.q_o),
        "q_4 evals:",
        abbreviate_vec(&composer.q_4),
        "q_c evals:",
        abbreviate_vec(&composer.q_c),
        "q_arith evals:",
        abbreviate_vec(&composer.q_arith),
        "q_range evals:",
        abbreviate_vec(&composer.q_range),
        "q_logic evals:",
        abbreviate_vec(&composer.q_logic),
        "q_fixed_group_add evals:",
        abbreviate_vec(&composer.q_fixed_group_add),
        "q_variable_group_add evals:",
        abbreviate_vec(
            &composer.q_variable_group_add
        ),
        "q_lookup evals:",
        abbreviate_vec(&composer.q_lookup),
        composer.public_inputs_sparse_store,
        "wires:",
        wires_debug_print::<F>(vec![
            &composer.w_l, &composer.w_r, &composer.w_o, &composer.w_4
        ]),
        "wires evals:",
        wires_values_debug_print(vec![
            &composer.w_l, &composer.w_r, &composer.w_o, &composer.w_4
        ], &composer.variables),
        "lookup table:",
        lookup_table_debug_print::<F>(&composer.lookup_table),
        "wire values:",
        var_map_debug_print(&composer.variables),
        composer.perm,
    );
    debug_str
}

/// pretty prints a lookup table for debugging
pub fn lookup_table_debug_print<F: PrimeField>(lt: &LookupTable<F>) -> String {
    fn abbreviate_field<F>(f: F) -> String
    where
        F: Field,
    {
        let whole_num_string = f
            .to_string()
            .split(|c| c == '(' || c == ')')
            .collect::<Vec<_>>()[1]
            .to_string();
        let len = whole_num_string.len();
        [&whole_num_string[..3], &whole_num_string[(len - 3)..]].join("##")
    }
    let debug_str = format!(
        "{}\n    {}",
        ["    column 1", "column 2", "column 3", "column 4"].join("      "),
        lt.0
            .iter()
            .map(|r| r
                .iter()
                .map(|f| abbreviate_field::<F>(*f))
                .collect::<Vec<_>>()
                .join("      "))
            .collect::<Vec<_>>()
            .join("\n    ")
    );

    debug_str
}