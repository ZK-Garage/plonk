// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.

//! Proof System Widgets

pub mod arithmetic;
pub mod ecc;
pub mod logic;
pub mod range;

use crate::proof_system::linearisation_poly::ProofEvaluations;
use crate::proof_system::permutation;
use crate::transcript::TranscriptProtocol;
use ark_ec::PairingEngine;
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::{univariate::DensePolynomial, Evaluations};
use ark_poly_commit::sonic_pc::SonicKZG10;
use ark_poly_commit::{PCCommitment, PolynomialCommitment};
use ark_serialize::*;

/// Gate Values
///
/// This data structures holds the wire values for a given gate.
#[derive(derivative::Derivative)]
#[derivative(Clone, Copy, Debug, Default, Eq, Hash, PartialEq)]
pub struct GateValues<F>
where
    F: Field,
{
    /// Left Value
    pub left: F,

    /// Right Value
    pub right: F,

    /// Output Value
    pub output: F,

    /// Fourth Value
    pub fourth: F,

    /// Next Left Value
    ///
    /// Only used in gates which use internal copy constraints.
    pub left_next: F,

    /// Next Right Value
    ///
    /// Only used in gates which use internal copy constraints.
    pub right_next: F,

    /// Next Fourth Value
    ///
    /// Only used in gates which use internal copy constraints.
    pub fourth_next: F,

    /// Left Wire Selector Weight
    pub left_selector: F,

    /// Right Wire Selector Weight
    pub right_selector: F,

    /// Constant Wire Selector Weight
    pub constant_selector: F,
}

/// Gate Constraint
pub trait GateConstraint<F>
where
    F: Field,
{
    /// Returns the coefficient of the quotient polynomial for this gate given
    /// an instantiation of the gate at `values` and a
    /// `separation_challenge` if this gate requires it for soundness.
    ///
    /// This method is an encoding of the polynomial constraint(s) that this
    /// gate represents whenever it is added to a circuit.
    fn constraints(separation_challenge: F, values: GateValues<F>) -> F {
        unimplemented!()
    }

    /// Computes the quotient polynomial term for the given gate type for the
    /// given value of `selector` instantiated with `separation_challenge` and
    /// `values`.
    fn quotient_term(
        selector: F,
        separation_challenge: F,
        values: GateValues<F>,
    ) -> F {
        selector * Self::constraints(separation_challenge, values)
    }

    /// Computes the linearisation polynomial term for the given gate type
    /// at the `selector_polynomial` instantiated with `separation_challenge`
    /// and `values`.
    fn linearisation_term(
        selector_polynomial: &DensePolynomial<F>,
        separation_challenge: F,
        values: GateValues<F>,
    ) -> DensePolynomial<F> {
        selector_polynomial * Self::constraints(separation_challenge, values)
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
        PC: PolynomialCommitment<F, DensePolynomial<F>>,
    {
        let coefficient = Self::constraints(
            separation_challenge,
            GateValues {
                left: evaluations.a_eval,
                right: evaluations.b_eval,
                output: evaluations.c_eval,
                fourth: evaluations.d_eval,
                left_next: evaluations.a_next_eval,
                right_next: evaluations.b_next_eval,
                fourth_next: evaluations.d_next_eval,
                left_selector: evaluations.q_l_eval,
                right_selector: evaluations.q_r_eval,
                constant_selector: evaluations.q_c_eval,
            },
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
    //Debug(bound = ""),
    //Eq(bound = ""),
    //PartialEq(bound = "")
)]
pub struct VerifierKey<F, PC>
where
    F: PrimeField,
    PC: PolynomialCommitment<F, DensePolynomial<F>>,
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
}

impl<F, PC> VerifierKey<F, PC>
where
    F: PrimeField,
    PC: PolynomialCommitment<F, DensePolynomial<F>>,
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
        q_arith: PC::Commitment,
        q_range: PC::Commitment,
        q_logic: PC::Commitment,
        q_fixed_group_add: PC::Commitment,
        q_variable_group_add: PC::Commitment,
        left_sigma: PC::Commitment,
        right_sigma: PC::Commitment,
        out_sigma: PC::Commitment,
        fourth_sigma: PC::Commitment,
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
    PC: PolynomialCommitment<F, DensePolynomial<F>>,
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
    Clone(bound = ""),
    Debug(bound = ""),
    Eq(bound = ""),
    PartialEq(bound = "")
)]
pub struct ProverKey<F>
where
    F: FftField,
{
    /// Circuit size
    pub(crate) n: usize,

    /// Arithmetic Prover Key
    pub(crate) arithmetic: arithmetic::ProverKey<F>,

    /// Range Gate Selector
    pub(crate) range_selector: (DensePolynomial<F>, Evaluations<F>),

    /// Logic Gate Selector
    pub(crate) logic_selector: (DensePolynomial<F>, Evaluations<F>),

    /// Fixed Group Addition Selector
    pub(crate) fixed_group_add_selector: (DensePolynomial<F>, Evaluations<F>),

    /// Variable Group Addition Selector
    pub(crate) variable_group_add_selector:
        (DensePolynomial<F>, Evaluations<F>),

    /// ProverKey for permutation checks
    pub(crate) permutation: permutation::ProverKey<F>,

    /// Pre-processes the 4n Evaluations for the vanishing polynomial, so
    /// they do not need to be computed at the proving stage.
    ///
    /// NOTE: With this, we can combine all parts of the quotient polynomial
    /// in their evaluation phase and divide by the quotient
    /// polynomial without having to perform IFFT
    pub(crate) v_h_coset_4n: Evaluations<F>,
}

impl<F> ProverKey<F>
where
    F: FftField,
{
    pub(crate) fn v_h_coset_4n(&self) -> &Evaluations<F> {
        &self.v_h_coset_4n
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
        q_arith: (DensePolynomial<F>, Evaluations<F>),
        q_range: (DensePolynomial<F>, Evaluations<F>),
        q_logic: (DensePolynomial<F>, Evaluations<F>),
        q_fixed_group_add: (DensePolynomial<F>, Evaluations<F>),
        q_variable_group_add: (DensePolynomial<F>, Evaluations<F>),
        left_sigma: (DensePolynomial<F>, Evaluations<F>),
        right_sigma: (DensePolynomial<F>, Evaluations<F>),
        out_sigma: (DensePolynomial<F>, Evaluations<F>),
        fourth_sigma: (DensePolynomial<F>, Evaluations<F>),
        linear_evaluations: Evaluations<F>,
        v_h_coset_4n: Evaluations<F>,
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
                q_arith,
            },
            range_selector: q_range,
            logic_selector: q_logic,
            fixed_group_add_selector: q_fixed_group_add,
            variable_group_add_selector: q_variable_group_add,
            permutation: permutation::ProverKey {
                left_sigma,
                right_sigma,
                out_sigma,
                fourth_sigma,
                linear_evaluations,
            },
            v_h_coset_4n,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_bls12_381::Fr as BlsScalar;
    use ark_bls12_381::G1Affine;
    use ark_ff::{Fp256, UniformRand};
    use ark_poly::polynomial::univariate::DensePolynomial;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, UVPolynomial};
    use ark_poly_commit::kzg10::Commitment;
    use rand_core::OsRng;

    fn rand_poly_eval(
        n: usize,
    ) -> (
        DensePolynomial<Fp256<ark_bls12_381::FrParameters>>,
        Evaluations<Fp256<ark_bls12_381::FrParameters>>,
    ) {
        let polynomial = DensePolynomial::rand(n, &mut OsRng);
        (polynomial, rand_evaluations(n))
    }

    fn rand_evaluations(
        n: usize,
    ) -> Evaluations<Fp256<ark_bls12_381::FrParameters>> {
        let domain: GeneralEvaluationDomain<
            Fp256<ark_bls12_381::FrParameters>,
        > = GeneralEvaluationDomain::new(4 * n).unwrap();
        let values: Vec<_> =
            (0..4 * n).map(|_| BlsScalar::rand(&mut OsRng)).collect();
        Evaluations::from_vec_and_domain(values, domain)
    }

    #[test]
    fn test_serialise_deserialise_prover_key() {
        let n = 1 << 11;

        let q_m = rand_poly_eval(n);
        let q_l = rand_poly_eval(n);
        let q_r = rand_poly_eval(n);
        let q_o = rand_poly_eval(n);
        let q_4 = rand_poly_eval(n);
        let q_c = rand_poly_eval(n);
        let q_arith = rand_poly_eval(n);
        let q_range = rand_poly_eval(n);
        let q_logic = rand_poly_eval(n);
        let q_fixed_group_add = rand_poly_eval(n);
        let q_variable_group_add = rand_poly_eval(n);

        let left_sigma = rand_poly_eval(n);
        let right_sigma = rand_poly_eval(n);
        let out_sigma = rand_poly_eval(n);
        let fourth_sigma = rand_poly_eval(n);

        let linear_evaluations = rand_evaluations(n);
        let v_h_coset_4n = rand_evaluations(n);

        let prover_key = ProverKey::from_polynomials_and_evals(
            n,
            q_m,
            q_l,
            q_r,
            q_o,
            q_4,
            q_c,
            q_arith,
            q_range,
            q_logic,
            q_fixed_group_add,
            q_variable_group_add,
            left_sigma,
            right_sigma,
            out_sigma,
            fourth_sigma,
            linear_evaluations,
            v_h_coset_4n,
        );

        let mut prover_key_bytes = vec![];
        prover_key
            .serialize_unchecked(&mut prover_key_bytes)
            .unwrap();

        let obtained_pk: ProverKey<Fp256<ark_bls12_381::FrParameters>> =
            ProverKey::deserialize_unchecked(prover_key_bytes.as_slice())
                .unwrap();

        assert!(prover_key == obtained_pk);
    }
    /*
        #[test]
        fn test_serialise_deserialise_verifier_key<E>()
        where
            E: PairingEngine,
        {
            let n = 2usize.pow(5);

            let q_m = Commitment::<E>(G1Affine::default());
            let q_l = Commitment(G1Affine::default());
            let q_r = Commitment(G1Affine::default());
            let q_o = Commitment(G1Affine::default());
            let q_4 = Commitment(G1Affine::default());
            let q_c = Commitment(G1Affine::default());
            let q_arith = Commitment(G1Affine::default());
            let q_range = Commitment(G1Affine::default());
            let q_logic = Commitment(G1Affine::default());
            let q_fixed_group_add = Commitment(G1Affine::default());
            let q_variable_group_add = Commitment(G1Affine::default());

            let left_sigma = Commitment(G1Affine::default());
            let right_sigma = Commitment(G1Affine::default());
            let out_sigma = Commitment(G1Affine::default());
            let fourth_sigma = Commitment(G1Affine::default());

            let verifier_key = VerifierKey::from_polynomial_commitments(
                n,
                q_m,
                q_l,
                q_r,
                q_o,
                q_4,
                q_c,
                q_arith,
                q_range,
                q_logic,
                q_fixed_group_add,
                q_variable_group_add,
                left_sigma,
                right_sigma,
                out_sigma,
                fourth_sigma,
            );

            let mut verifier_key_bytes = vec![];
            verifier_key
                .serialize_unchecked(&mut verifier_key_bytes)
                .unwrap();

            type Fr = <Bls12_381 as PairingEngine>::Fr;
            let obtained_vk: VerifierKey<
                Fr,
                SonicKZG10<Bls12_381, DensePolynomial<Fr>>,
            > = VerifierKey::deserialize_unchecked(verifier_key_bytes.as_slice())
                .unwrap();

            //assert!(verifier_key == obtained_vk);
        }
    */
}
