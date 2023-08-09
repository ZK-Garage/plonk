// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) ZK-Garage. All rights reserved.

//! Proof System Widgets

pub mod arithmetic;
pub mod ecc;
pub mod logic;
pub mod lookup;
pub mod range;

use crate::{
    commitment::HomomorphicCommitment,
    lookup::MultiSet,
    proof_system::{
        linearisation_poly::CustomEvaluations,
        linearisation_poly::ProofEvaluations, permutation,
    },
    transcript::TranscriptProtocol,
};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, Evaluations};
use ark_serialize::*;

/// Set of values needed for a custom gate
pub trait CustomValues<F>
where
    F: PrimeField,
{
    /// Constructs gate-specific values struct from the set of evaluations
    /// `CustomEvaluations`
    fn from_evaluations(custom_evals: &CustomEvaluations<F>) -> Self;
}

/// Witness Values
///
/// This data structures holds the wire values for a given gate.
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct WitnessValues<F>
where
    F: PrimeField,
{
    /// Left Value
    pub a_val: F,

    /// Right Value
    pub b_val: F,

    /// Output Value
    pub c_val: F,

    /// Fourth Value
    pub d_val: F,
}

/// Gate Constraint
pub trait GateConstraint<F>
where
    F: PrimeField,
{
    /// Custom values needed for the gate
    type CustomVals: CustomValues<F>;

    /// Returns the coefficient of the quotient polynomial for this gate given
    /// an instantiation of the gate at `values` and a
    /// `separation_challenge` if this gate requires it for soundness.
    ///
    /// This method is an encoding of the polynomial constraint(s) that this
    /// gate represents whenever it is added to a circuit.
    fn constraints(
        separation_challenge: F,
        wit_vals: WitnessValues<F>,
        custom_vals: Self::CustomVals,
    ) -> F;

    /// Computes the quotient polynomial term for the given gate type for the
    /// given value of `selector` instantiated with `separation_challenge` and
    /// `values`.
    fn quotient_term(
        selector: F,
        separation_challenge: F,
        wit_vals: WitnessValues<F>,
        custom_vals: Self::CustomVals,
    ) -> F {
        selector
            * Self::constraints(separation_challenge, wit_vals, custom_vals)
    }

    /// Computes the linearisation polynomial term for the given gate type
    /// at the `selector_polynomial` instantiated with `separation_challenge`
    /// and `values`.
    fn linearisation_term(
        selector_polynomial: &DensePolynomial<F>,
        separation_challenge: F,
        wit_vals: WitnessValues<F>,
        custom_vals: Self::CustomVals,
    ) -> DensePolynomial<F> {
        selector_polynomial
            * Self::constraints(separation_challenge, wit_vals, custom_vals)
    }

    /// Extends `scalars` and `points` to build the linearisation commitment
    /// with the given instantiation of `evaluations` and
    /// `separation_challenge`.
    fn extend_linearisation_commitment<PC>(
        selector_commitment: &PC::Commitment,
        separation_challenge: F,
        evaluations: &ProofEvaluations<F>,
        scalars: &mut Vec<F>,
        points: &mut Vec<PC::Commitment>,
    ) where
        PC: HomomorphicCommitment<F>,
    {
        let coefficient = Self::constraints(
            separation_challenge,
            WitnessValues {
                a_val: evaluations.wire_evals.a_eval,
                b_val: evaluations.wire_evals.b_eval,
                c_val: evaluations.wire_evals.c_eval,
                d_val: evaluations.wire_evals.d_eval,
            },
            Self::CustomVals::from_evaluations(&evaluations.custom_evals),
        );
        scalars.push(coefficient);
        points.push(selector_commitment.clone());
    }
}

/// PLONK circuit Verification Key.
///
/// This structure is used by the Verifier in order to verify a
/// [`Proof`](super::Proof).
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(
        bound = "arithmetic::VerifierKey<F,PC>: core::fmt::Debug, PC::Commitment: core::fmt::Debug"
    ),
    Eq(bound = "arithmetic::VerifierKey<F,PC>: Eq, PC::Commitment: Eq"),
    PartialEq(
        bound = "arithmetic::VerifierKey<F,PC>: PartialEq, PC::Commitment: PartialEq"
    )
)]
pub struct VerifierKey<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    /// Circuit size (not padded to a power of two).
    pub(crate) n: usize,

    /// Arithmetic Verifier Key
    pub(crate) arithmetic: arithmetic::VerifierKey<F, PC>,

    /// Range Gate Selector Commitment
    pub(crate) range_selector_commitment: PC::Commitment,

    /// Logic Gate Selector Commitment
    pub(crate) logic_selector_commitment: PC::Commitment,

    /// Fixed Group Addition Selector Commitment
    pub(crate) fixed_group_add_selector_commitment: PC::Commitment,

    /// Variable Group Addition Selector Commitment
    pub(crate) variable_group_add_selector_commitment: PC::Commitment,

    /// VerifierKey for permutation checks
    pub(crate) permutation: permutation::VerifierKey<PC::Commitment>,

    /// VerifierKey for Lookup Gate
    pub(crate) lookup: lookup::VerifierKey<F, PC>,
}

impl<F, PC> VerifierKey<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    /// Constructs a [`VerifierKey`] from the widget VerifierKey's that are
    /// constructed based on the selector polynomial commitments and the
    /// sigma polynomial commitments.
    pub(crate) fn from_polynomial_commitments(
        n: usize,
        q_m: PC::Commitment,
        q_l: PC::Commitment,
        q_r: PC::Commitment,
        q_o: PC::Commitment,
        q_4: PC::Commitment,
        q_c: PC::Commitment,
        q_hl: PC::Commitment,
        q_hr: PC::Commitment,
        q_h4: PC::Commitment,
        q_arith: PC::Commitment,
        q_range: PC::Commitment,
        q_logic: PC::Commitment,
        q_lookup: PC::Commitment,
        q_fixed_group_add: PC::Commitment,
        q_variable_group_add: PC::Commitment,
        left_sigma: PC::Commitment,
        right_sigma: PC::Commitment,
        out_sigma: PC::Commitment,
        fourth_sigma: PC::Commitment,
        table_1: PC::Commitment,
        table_2: PC::Commitment,
        table_3: PC::Commitment,
        table_4: PC::Commitment,
    ) -> Self {
        Self {
            n,
            arithmetic: arithmetic::VerifierKey {
                q_m,
                q_l,
                q_r,
                q_o,
                q_4,
                q_c,
                q_hl,
                q_hr,
                q_h4,
                q_arith,
            },
            range_selector_commitment: q_range,
            logic_selector_commitment: q_logic,
            fixed_group_add_selector_commitment: q_fixed_group_add,
            variable_group_add_selector_commitment: q_variable_group_add,
            permutation: permutation::VerifierKey {
                left_sigma,
                right_sigma,
                out_sigma,
                fourth_sigma,
            },
            lookup: lookup::VerifierKey {
                q_lookup,
                table_1,
                table_2,
                table_3,
                table_4,
            },
        }
    }

    /// Returns the Circuit size padded to the next power of two.
    pub fn padded_circuit_size(&self) -> usize {
        self.n.next_power_of_two()
    }
}

impl<F, PC> VerifierKey<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    /// Adds the circuit description to the transcript.
    pub(crate) fn seed_transcript<T>(&self, transcript: &mut T)
    where
        T: TranscriptProtocol,
    {
        transcript.append(b"q_m", &self.arithmetic.q_m);
        transcript.append(b"q_l", &self.arithmetic.q_l);
        transcript.append(b"q_r", &self.arithmetic.q_r);
        transcript.append(b"q_o", &self.arithmetic.q_o);
        transcript.append(b"q_c", &self.arithmetic.q_c);
        transcript.append(b"q_4", &self.arithmetic.q_4);
        transcript.append(b"q_hl", &self.arithmetic.q_hl);
        transcript.append(b"q_hr", &self.arithmetic.q_hr);
        transcript.append(b"q_h4", &self.arithmetic.q_h4);
        transcript.append(b"q_arith", &self.arithmetic.q_arith);
        transcript.append(b"q_range", &self.range_selector_commitment);
        transcript.append(b"q_logic", &self.logic_selector_commitment);
        transcript.append(
            b"q_variable_group_add",
            &self.variable_group_add_selector_commitment,
        );
        transcript.append(
            b"q_fixed_group_add",
            &self.fixed_group_add_selector_commitment,
        );
        transcript.append(b"left_sigma", &self.permutation.left_sigma);
        transcript.append(b"right_sigma", &self.permutation.right_sigma);
        transcript.append(b"out_sigma", &self.permutation.out_sigma);
        transcript.append(b"fourth_sigma", &self.permutation.fourth_sigma);
        transcript.circuit_domain_sep(self.n as u64);
    }
}

/// PLONK circuit Proving Key.
///
/// This structure is used by the Prover in order to construct a
/// [`Proof`](crate::proof_system::Proof).
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone(bound = "lookup::ProverKey<F>: Clone"),
    Debug(bound = "lookup::ProverKey<F>: std::fmt::Debug"),
    Eq(bound = "lookup::ProverKey<F>: Eq"),
    PartialEq(bound = "lookup::ProverKey<F>: PartialEq")
)]
pub struct ProverKey<F>
where
    F: PrimeField,
{
    /// Circuit size
    pub(crate) n: usize,

    /// Arithmetic Prover Key
    pub(crate) arithmetic: arithmetic::ProverKey<F>,

    /// Range Gate Selector
    pub(crate) range_selector: (DensePolynomial<F>, Evaluations<F>),

    /// Logic Gate Selector
    pub(crate) logic_selector: (DensePolynomial<F>, Evaluations<F>),

    /// Lookup selector
    pub(crate) lookup: lookup::ProverKey<F>,

    /// Fixed Group Addition Selector
    pub(crate) fixed_group_add_selector: (DensePolynomial<F>, Evaluations<F>),

    /// Variable Group Addition Selector
    pub(crate) variable_group_add_selector:
        (DensePolynomial<F>, Evaluations<F>),

    /// ProverKey for permutation checks
    pub(crate) permutation: permutation::ProverKey<F>,

    /// Pre-processes the 8n Evaluations for the vanishing polynomial, so
    /// they do not need to be computed at the proving stage.
    ///
    /// NOTE: With this, we can combine all parts of the quotient polynomial
    /// in their evaluation phase and divide by the quotient
    /// polynomial without having to perform IFFT
    pub(crate) v_h_coset_8n: Evaluations<F>,
}

impl<F> ProverKey<F>
where
    F: PrimeField,
{
    pub(crate) fn v_h_coset_8n(&self) -> &Evaluations<F> {
        &self.v_h_coset_8n
    }

    /// Constructs a [`ProverKey`] from the widget ProverKey's that are
    /// constructed based on the selector polynomials and the
    /// sigma polynomials and it's evaluations.
    pub(crate) fn from_polynomials_and_evals(
        n: usize,
        q_m: (DensePolynomial<F>, Evaluations<F>),
        q_l: (DensePolynomial<F>, Evaluations<F>),
        q_r: (DensePolynomial<F>, Evaluations<F>),
        q_o: (DensePolynomial<F>, Evaluations<F>),
        q_4: (DensePolynomial<F>, Evaluations<F>),
        q_c: (DensePolynomial<F>, Evaluations<F>),
        q_hl: (DensePolynomial<F>, Evaluations<F>),
        q_hr: (DensePolynomial<F>, Evaluations<F>),
        q_h4: (DensePolynomial<F>, Evaluations<F>),
        q_arith: (DensePolynomial<F>, Evaluations<F>),
        q_range: (DensePolynomial<F>, Evaluations<F>),
        q_logic: (DensePolynomial<F>, Evaluations<F>),
        q_lookup: (DensePolynomial<F>, Evaluations<F>),
        q_fixed_group_add: (DensePolynomial<F>, Evaluations<F>),
        q_variable_group_add: (DensePolynomial<F>, Evaluations<F>),
        left_sigma: (DensePolynomial<F>, Evaluations<F>),
        right_sigma: (DensePolynomial<F>, Evaluations<F>),
        out_sigma: (DensePolynomial<F>, Evaluations<F>),
        fourth_sigma: (DensePolynomial<F>, Evaluations<F>),
        linear_evaluations: Evaluations<F>,
        v_h_coset_8n: Evaluations<F>,
        table_1: MultiSet<F>,
        table_2: MultiSet<F>,
        table_3: MultiSet<F>,
        table_4: MultiSet<F>,
    ) -> Self {
        Self {
            n,
            arithmetic: arithmetic::ProverKey {
                q_m,
                q_l,
                q_r,
                q_o,
                q_4,
                q_c,
                q_hl,
                q_hr,
                q_h4,
                q_arith,
            },
            range_selector: q_range,
            logic_selector: q_logic,
            fixed_group_add_selector: q_fixed_group_add,
            variable_group_add_selector: q_variable_group_add,
            lookup: lookup::ProverKey {
                q_lookup,
                table_1,
                table_2,
                table_3,
                table_4,
            },
            permutation: permutation::ProverKey {
                left_sigma,
                right_sigma,
                out_sigma,
                fourth_sigma,
                linear_evaluations,
            },
            v_h_coset_8n,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::batch_test;
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_ec::models::TEModelParameters;
    use ark_poly::polynomial::univariate::DensePolynomial;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, UVPolynomial};
    use rand_core::OsRng;

    fn rand_poly_eval<F>(n: usize) -> (DensePolynomial<F>, Evaluations<F>)
    where
        F: PrimeField,
    {
        let polynomial = DensePolynomial::rand(n, &mut OsRng);
        (polynomial, rand_evaluations(n))
    }

    fn rand_evaluations<F>(n: usize) -> Evaluations<F>
    where
        F: PrimeField,
    {
        let domain = GeneralEvaluationDomain::new(4 * n).unwrap();
        let values: Vec<_> = (0..8 * n).map(|_| F::rand(&mut OsRng)).collect();
        Evaluations::from_vec_and_domain(values, domain)
    }

    fn rand_multiset<F>(n: usize) -> MultiSet<F>
    where
        F: PrimeField,
    {
        let mut rng = OsRng;
        core::iter::from_fn(|| Some(F::rand(&mut rng)))
            .take(n)
            .collect()
    }

    #[test]
    fn test_serialise_deserialise_prover_key() {
        type F = ark_bls12_381::Fr;
        let n = 1 << 11;

        let q_m = rand_poly_eval(n);
        let q_l = rand_poly_eval(n);
        let q_r = rand_poly_eval(n);
        let q_o = rand_poly_eval(n);
        let q_4 = rand_poly_eval(n);
        let q_c = rand_poly_eval(n);
        let q_hl = rand_poly_eval(n);
        let q_hr = rand_poly_eval(n);
        let q_h4 = rand_poly_eval(n);
        let q_arith = rand_poly_eval(n);
        let q_range = rand_poly_eval(n);
        let q_logic = rand_poly_eval(n);
        let q_lookup = rand_poly_eval(n);
        let q_fixed_group_add = rand_poly_eval(n);
        let q_variable_group_add = rand_poly_eval(n);

        let left_sigma = rand_poly_eval(n);
        let right_sigma = rand_poly_eval(n);
        let out_sigma = rand_poly_eval(n);
        let fourth_sigma = rand_poly_eval(n);

        let linear_evaluations = rand_evaluations(n);
        let v_h_coset_8n = rand_evaluations(n);
        let table_1 = rand_multiset(n);
        let table_2 = rand_multiset(n);
        let table_3 = rand_multiset(n);
        let table_4 = rand_multiset(n);

        let prover_key = ProverKey::from_polynomials_and_evals(
            n,
            q_m,
            q_l,
            q_r,
            q_o,
            q_4,
            q_c,
            q_hl,
            q_hr,
            q_h4,
            q_arith,
            q_range,
            q_logic,
            q_lookup,
            q_fixed_group_add,
            q_variable_group_add,
            left_sigma,
            right_sigma,
            out_sigma,
            fourth_sigma,
            linear_evaluations,
            v_h_coset_8n,
            table_1,
            table_2,
            table_3,
            table_4,
        );

        let mut prover_key_bytes = vec![];
        prover_key
            .serialize_unchecked(&mut prover_key_bytes)
            .unwrap();

        let obtained_pk: ProverKey<F> =
            ProverKey::deserialize_unchecked(prover_key_bytes.as_slice())
                .unwrap();

        assert_eq!(prover_key, obtained_pk);
    }

    #[allow(clippy::extra_unused_type_parameters)]
    fn test_serialise_deserialise_verifier_key<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
        VerifierKey<F, PC>: PartialEq,
    {
        let n = 2usize.pow(5);

        let q_m = PC::Commitment::default();
        let q_l = PC::Commitment::default();
        let q_r = PC::Commitment::default();
        let q_o = PC::Commitment::default();
        let q_4 = PC::Commitment::default();
        let q_c = PC::Commitment::default();

        let q_hl = PC::Commitment::default();
        let q_hr = PC::Commitment::default();
        let q_h4 = PC::Commitment::default();

        let q_arith = PC::Commitment::default();
        let q_range = PC::Commitment::default();
        let q_logic = PC::Commitment::default();
        let q_lookup = PC::Commitment::default();
        let q_fixed_group_add = PC::Commitment::default();
        let q_variable_group_add = PC::Commitment::default();

        let left_sigma = PC::Commitment::default();
        let right_sigma = PC::Commitment::default();
        let out_sigma = PC::Commitment::default();
        let fourth_sigma = PC::Commitment::default();

        let table_1 = PC::Commitment::default();
        let table_2 = PC::Commitment::default();
        let table_3 = PC::Commitment::default();
        let table_4 = PC::Commitment::default();

        let verifier_key = VerifierKey::<F, PC>::from_polynomial_commitments(
            n,
            q_m,
            q_l,
            q_r,
            q_o,
            q_4,
            q_c,
            q_hl,
            q_hr,
            q_h4,
            q_arith,
            q_range,
            q_logic,
            q_lookup,
            q_fixed_group_add,
            q_variable_group_add,
            left_sigma,
            right_sigma,
            out_sigma,
            fourth_sigma,
            table_1,
            table_2,
            table_3,
            table_4,
        );

        let mut verifier_key_bytes = vec![];
        verifier_key
            .serialize_unchecked(&mut verifier_key_bytes)
            .unwrap();

        let obtained_vk: VerifierKey<F, PC> =
            VerifierKey::deserialize_unchecked(verifier_key_bytes.as_slice())
                .unwrap();

        assert!(verifier_key == obtained_vk);
    }

    // Test for Bls12_381
    batch_test!(
        [test_serialise_deserialise_verifier_key],
        [] => (
            Bls12_381, ark_ed_on_bls12_381::EdwardsParameters      )
    );

    // Test for Bls12_377
    batch_test!(
        [test_serialise_deserialise_verifier_key],
        [] => (
            Bls12_377, ark_ed_on_bls12_377::EdwardsParameters       )
    );
}
