// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) DUSK NETWORK. All rights reserved.
/// XXX: Doc this

/// XXX: Doc this
pub mod arithmetic;
/// XXX: Doc this
pub mod ecc;
/// XXX: Doc this
pub mod logic;
/// XXX: Doc this
pub mod permutation;
/// XXX: Doc this
pub mod range;

use crate::transcript::TranscriptProtocol;
use ark_ec::{PairingEngine, TEModelParameters};
use ark_ff::PrimeField;
use ark_poly::{univariate::DensePolynomial, Evaluations};
use ark_poly_commit::sonic_pc::Commitment;
use ark_serialize::*;

/// PLONK circuit Verification Key.
///
/// This structure is used by the Verifier in order to verify a
/// [`Proof`](super::Proof).
#[derive(
    Debug, PartialEq, Eq, Copy, Clone, CanonicalDeserialize, CanonicalSerialize,
)]
pub struct VerifierKey<
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
> {
    /// Circuit size (not padded to a power of two).
    pub(crate) n: usize,
    /// VerifierKey for arithmetic gates
    pub(crate) arithmetic: arithmetic::VerifierKey<E>,
    /// VerifierKey for logic gates
    pub(crate) logic: logic::VerifierKey<E>,
    /// VerifierKey for range gates
    pub(crate) range: range::VerifierKey<E>,
    /// VerifierKey for fixed base curve addition gates
    pub(crate) fixed_base: ecc::scalar_mul::fixed_base::VerifierKey<E, P>,
    /// VerifierKey for variable base curve addition gates
    pub(crate) variable_base: ecc::curve_addition::VerifierKey<E, P>,
    /// VerifierKey for permutation checks
    pub(crate) permutation: permutation::VerifierKey<E>,
}

impl<E: PairingEngine, P: TEModelParameters<BaseField = E::Fr>>
    VerifierKey<E, P>
{
    /// Returns the Circuit size padded to the next power of two.
    pub fn padded_circuit_size(&self) -> usize {
        self.n.next_power_of_two()
    }

    /// Constructs a [`VerifierKey`] from the widget VerifierKey's that are
    /// constructed based on the selector polynomial commitments and the
    /// sigma polynomial commitments.
    pub(crate) fn from_polynomial_commitments(
        n: usize,
        q_m: Commitment<E>,
        q_l: Commitment<E>,
        q_r: Commitment<E>,
        q_o: Commitment<E>,
        q_4: Commitment<E>,
        q_c: Commitment<E>,
        q_arith: Commitment<E>,
        q_logic: Commitment<E>,
        q_range: Commitment<E>,
        q_fixed_group_add: Commitment<E>,
        q_variable_group_add: Commitment<E>,
        left_sigma: Commitment<E>,
        right_sigma: Commitment<E>,
        out_sigma: Commitment<E>,
        fourth_sigma: Commitment<E>,
    ) -> VerifierKey<E, P> {
        let arithmetic = arithmetic::VerifierKey {
            q_m,
            q_l,
            q_r,
            q_o,
            q_4,
            q_c,
            q_arith,
        };
        let logic = logic::VerifierKey { q_c, q_logic };
        let range = range::VerifierKey { q_range };
        let fixed_base = ecc::scalar_mul::fixed_base::VerifierKey::new(
            q_l,
            q_r,
            q_fixed_group_add,
        );

        let variable_base =
            ecc::curve_addition::VerifierKey::new(q_variable_group_add);

        let permutation = permutation::VerifierKey::new(
            left_sigma,
            right_sigma,
            out_sigma,
            fourth_sigma,
        );

        VerifierKey {
            n,
            arithmetic,
            logic,
            range,
            fixed_base,
            variable_base,
            permutation,
        }
    }
}

impl<E: PairingEngine, P: TEModelParameters<BaseField = E::Fr>>
    VerifierKey<E, P>
{
    /// Adds the circuit description to the transcript
    pub(crate) fn seed_transcript<T: TranscriptProtocol<E>>(
        &self,
        transcript: &mut T,
    ) {
        transcript.append_commitment(b"q_m", &self.arithmetic.q_m);
        transcript.append_commitment(b"q_l", &self.arithmetic.q_l);
        transcript.append_commitment(b"q_r", &self.arithmetic.q_r);
        transcript.append_commitment(b"q_o", &self.arithmetic.q_o);
        transcript.append_commitment(b"q_c", &self.arithmetic.q_c);
        transcript.append_commitment(b"q_4", &self.arithmetic.q_4);
        transcript.append_commitment(b"q_arith", &self.arithmetic.q_arith);
        transcript.append_commitment(b"q_range", &self.range.q_range);
        transcript.append_commitment(b"q_logic", &self.logic.q_logic);
        transcript.append_commitment(
            b"q_variable_group_add",
            &self.variable_base.q_variable_group_add,
        );
        transcript.append_commitment(
            b"q_fixed_group_add",
            &self.fixed_base.q_fixed_group_add,
        );

        transcript
            .append_commitment(b"left_sigma", &self.permutation.left_sigma);
        transcript
            .append_commitment(b"right_sigma", &self.permutation.right_sigma);
        transcript.append_commitment(b"out_sigma", &self.permutation.out_sigma);
        transcript
            .append_commitment(b"fourth_sigma", &self.permutation.fourth_sigma);

        // Append circuit size to transcript
        transcript.circuit_domain_sep(self.n as u64);
    }
}

/// PLONK circuit Proving Key.
///
/// This structure is used by the Prover in order to construct a
/// [`Proof`](crate::proof_system::Proof).
#[derive(
    Debug, PartialEq, Eq, Clone, CanonicalSerialize, CanonicalDeserialize,
)]
pub struct ProverKey<F: PrimeField, P: TEModelParameters<BaseField = F>> {
    /// Circuit size
    pub(crate) n: usize,
    /// ProverKey for arithmetic gate
    pub(crate) arithmetic: arithmetic::ProverKey<F>,
    /// ProverKey for logic gate
    pub(crate) logic: logic::ProverKey<F>,
    /// ProverKey for range gate
    pub(crate) range: range::ProverKey<F>,
    /// ProverKey for fixed base curve addition gates
    pub(crate) fixed_base: ecc::scalar_mul::fixed_base::ProverKey<F, P>,
    /// ProverKey for variable base curve addition gates
    pub(crate) variable_base: ecc::curve_addition::ProverKey<F, P>,
    /// ProverKey for permutation checks
    pub(crate) permutation: permutation::ProverKey<F>,
    // Pre-processes the 4n Evaluations for the vanishing polynomial, so
    // they do not need to be computed at the proving stage.
    // Note: With this, we can combine all parts of the quotient polynomial
    // in their evaluation phase and divide by the quotient
    // polynomial without having to perform IFFT
    pub(crate) v_h_coset_4n: Evaluations<F>,
}

impl<F: PrimeField, P: TEModelParameters<BaseField = F>> ProverKey<F, P> {
    /// Returns the number of [`Polynomial`]s contained in a ProverKey.
    fn num_polys() -> usize {
        15
    }

    /// Returns the number of [`Evaluations`] contained in a ProverKey.
    fn num_evals() -> usize {
        17
    }

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
        q_logic: (DensePolynomial<F>, Evaluations<F>),
        q_range: (DensePolynomial<F>, Evaluations<F>),
        q_fixed_group_add: (DensePolynomial<F>, Evaluations<F>),
        q_variable_group_add: (DensePolynomial<F>, Evaluations<F>),
        left_sigma: (DensePolynomial<F>, Evaluations<F>),
        right_sigma: (DensePolynomial<F>, Evaluations<F>),
        out_sigma: (DensePolynomial<F>, Evaluations<F>),
        fourth_sigma: (DensePolynomial<F>, Evaluations<F>),
        linear_evaluations: Evaluations<F>,
        v_h_coset_4n: Evaluations<F>,
    ) -> ProverKey<F, P> {
        let arithmetic = arithmetic::ProverKey {
            q_m,
            q_l: q_l.clone(),
            q_r: q_r.clone(),
            q_o,
            q_4,
            q_c: q_c.clone(),
            q_arith,
        };
        let logic = logic::ProverKey {
            q_c: q_c.clone(),
            q_logic,
        };
        let range = range::ProverKey { q_range };
        let fixed_base = ecc::scalar_mul::fixed_base::ProverKey::new(
            q_l,
            q_r,
            q_c,
            q_fixed_group_add,
        );

        let variable_base = ecc::curve_addition::ProverKey::new(
            q_variable_group_add.0,
            q_variable_group_add.1,
        );

        let permutation = permutation::ProverKey {
            left_sigma,
            right_sigma,
            out_sigma,
            fourth_sigma,
            linear_evaluations,
        };

        ProverKey {
            n,
            arithmetic,
            logic,
            range,
            fixed_base,
            variable_base,
            permutation,
            v_h_coset_4n,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use ark_bls12_381::Fr as BlsScalar;
    use ark_ff::UniformRand;
    use ark_poly::polynomial::univariate::DensePolynomial;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain, UVPolynomial};
    use rand_core::OsRng;

    fn rand_poly_eval<F: PrimeField>(
        n: usize,
    ) -> (DensePolynomial<F>, Evaluations<F>) {
        let polynomial = DensePolynomial::rand(n, &mut OsRng);
        (polynomial, rand_evaluations(n))
    }

    fn rand_evaluations<F: PrimeField>(n: usize) -> Evaluations<F> {
        let domain = GeneralEvaluationDomain::new(4 * n).unwrap();
        let values: Vec<_> =
            (0..4 * n).map(|_| BlsScalar::rand(&mut OsRng)).collect();
        let evaluations = Evaluations::from_vec_and_domain(values, domain);
        evaluations
    }

    #[test]
    fn test_serialise_deserialise_prover_key() {
        let n = 1 << 11;

        let q_m = rand_poly_eval(n);
        let q_l = rand_poly_eval(n);
        let q_r = rand_poly_eval(n);
        let q_o = rand_poly_eval(n);
        let q_c = rand_poly_eval(n);
        let q_4 = rand_poly_eval(n);
        let q_arith = rand_poly_eval(n);

        let q_logic = rand_poly_eval(n);

        let q_range = rand_poly_eval(n);

        let q_fixed_group_add = rand_poly_eval(n);

        let q_variable_group_add = rand_poly_eval(n);

        let left_sigma = rand_poly_eval(n);
        let right_sigma = rand_poly_eval(n);
        let out_sigma = rand_poly_eval(n);
        let fourth_sigma = rand_poly_eval(n);
        let linear_evaluations = rand_evaluations(n);

        let v_h_coset_4n = rand_evaluations(n);

        let arithmetic = arithmetic::ProverKey {
            q_m,
            q_l: q_l.clone(),
            q_r: q_r.clone(),
            q_o,
            q_c: q_c.clone(),
            q_4,
            q_arith,
        };

        let logic = logic::ProverKey {
            q_logic,
            q_c: q_c.clone(),
        };

        let range = range::ProverKey { q_range };

        let fixed_base = ecc::scalar_mul::fixed_base::ProverKey {
            q_fixed_group_add,
            q_l,
            q_r,
            q_c,
        };

        let permutation = permutation::ProverKey {
            left_sigma,
            right_sigma,
            out_sigma,
            fourth_sigma,
            linear_evaluations,
        };

        let variable_base = ecc::curve_addition::ProverKey {
            q_variable_group_add,
        };

        let prover_key = ProverKey {
            n,
            arithmetic,
            logic,
            fixed_base,
            range,
            variable_base,
            permutation,
            v_h_coset_4n,
        };

        let prover_key_bytes = prover_key.to_var_bytes();
        let pk = ProverKey::from_slice(&prover_key_bytes).unwrap();

        assert_eq!(pk, prover_key);
        assert_eq!(pk.to_var_bytes(), prover_key.to_var_bytes());
    }

    #[test]
    fn test_serialise_deserialise_verifier_key() {
        use ark_poly_commit::Commitment;
        use dusk_bls12_381::G1Affine;

        let n = 2usize.pow(5);

        let q_m = Commitment(G1Affine::generator());
        let q_l = Commitment(G1Affine::generator());
        let q_r = Commitment(G1Affine::generator());
        let q_o = Commitment(G1Affine::generator());
        let q_c = Commitment(G1Affine::generator());
        let q_4 = Commitment(G1Affine::generator());
        let q_arith = Commitment(G1Affine::generator());

        let q_range = Commitment(G1Affine::generator());

        let q_fixed_group_add = Commitment(G1Affine::generator());
        let q_variable_group_add = Commitment(G1Affine::generator());

        let q_logic = Commitment(G1Affine::generator());

        let left_sigma = Commitment(G1Affine::generator());
        let right_sigma = Commitment(G1Affine::generator());
        let out_sigma = Commitment(G1Affine::generator());
        let fourth_sigma = Commitment(G1Affine::generator());

        let arithmetic = arithmetic::VerifierKey {
            q_m,
            q_l,
            q_r,
            q_o,
            q_c,
            q_4,
            q_arith,
        };

        let logic = logic::VerifierKey { q_logic, q_c };

        let range = range::VerifierKey { q_range };

        let fixed_base = ecc::scalar_mul::fixed_base::VerifierKey {
            q_fixed_group_add,
            q_l,
            q_r,
        };
        let variable_base = ecc::curve_addition::VerifierKey {
            q_variable_group_add,
        };

        let permutation = permutation::VerifierKey {
            left_sigma,
            right_sigma,
            out_sigma,
            fourth_sigma,
        };

        let verifier_key = VerifierKey {
            n,
            arithmetic,
            logic,
            range,
            fixed_base,
            variable_base,
            permutation,
        };

        let verifier_key_bytes = verifier_key.to_bytes();
        let got = VerifierKey::from_bytes(&verifier_key_bytes).unwrap();

        assert_eq!(got, verifier_key);
    }
}
