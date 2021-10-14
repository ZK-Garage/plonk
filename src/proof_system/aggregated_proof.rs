use ark_ec::PairingEngine;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly::{Evaluations, GeneralEvaluationDomain};
use ark_poly_commit::kzg10::Commitment;
use merlin::Transcript;

/// Proof that multiple polynomials were correctly evaluated at a point `z`,
/// each producing their respective evaluated points p_i(z).
#[derive(Debug)]
pub(crate) struct PCAggregateProof<E: PairingEngine> {
    /// This is a commitment to the aggregated witness polynomial.
    pub(crate) commitment_to_witness: Commitment<E>,
    /// These are the results of the evaluating each polynomial at the
    /// point `z`.
    pub(crate) evaluated_points: Vec<E::Fr>,
    /// These are the commitments to the polynomials which you want to
    /// prove a statement about.
    pub(crate) commitments_to_polynomials: Vec<Commitment<E>>,
}

impl<E: PairingEngine> PCAggregateProof<E> {
    /// Initialises an `AggregatedProof` with the commitment to the witness.
    pub(crate) fn with_witness(witness: Commitment<E>) -> PCAggregateProof<E> {
        PCAggregateProof {
            commitment_to_witness: witness,
            evaluated_points: Vec::new(),
            commitments_to_polynomials: Vec::new(),
        }
    }

    /// Adds an evaluated point with the commitment to the polynomial which
    /// produced it.
    pub(crate) fn add_part(&mut self, part: (E::Fr, Commitment<E>)) {
        self.evaluated_points.push(part.0);
        self.commitments_to_polynomials.push(part.1);
    }

    /// Flattens an `PCAggregateProof` into a `Proof`.
    /// The transcript must have the same view as the transcript that was
    /// used to aggregate the witness in the proving stage.
    pub(crate) fn flatten(&self, transcript: &mut Transcript) -> PCProof {
        let challenge = transcript.challenge_scalar(b"aggregate_witness");
        let powers = crate::util::powers_of(
            &challenge,
            self.commitments_to_polynomials.len() - 1,
        );

        let flattened_poly_commitments_iter = self
            .commitments_to_polynomials
            .par_iter()
            .zip(powers.par_iter());
        let flattened_poly_evaluations_iter =
            self.evaluated_points.par_iter().zip(powers.par_iter());

        // Flattened polynomial commitments using challenge
        let flattened_poly_commitments: G1Projective =
            flattened_poly_commitments_iter
                .map(|(poly, challenge)| poly.0 * challenge)
                .sum();
        // Flattened evaluation points
        let flattened_poly_evaluations: BlsScalar =
            flattened_poly_evaluations_iter
                .map(|(eval, challenge)| eval * challenge)
                .sum();

        Proof {
            commitment_to_witness: self.commitment_to_witness,
            evaluated_point: flattened_poly_evaluations,
            commitment_to_polynomial: Commitment::from(
                flattened_poly_commitments,
            ),
        }
    }
}
