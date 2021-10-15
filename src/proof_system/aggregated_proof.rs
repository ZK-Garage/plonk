use ark_ec::{AffineCurve, PairingEngine, TEModelParameters};
use ark_poly_commit::kzg10::{Commitment, Proof as KZGProof};
use core::marker::PhantomData;
use merlin::Transcript;

use crate::transcript::TranscriptProtocol;

/// Proof that multiple polynomials were correctly evaluated at a point `z`,
/// each producing their respective evaluated points p_i(z).
#[derive(Debug)]
pub(crate) struct PCAggregateProof<
    E: PairingEngine<Fr = P::BaseField>,
    P: TEModelParameters,
> {
    /// This is a commitment to the aggregated witness polynomial.
    /// The aggregate witness polynomial is a linear combination of the
    /// witness polynomials (p_i(X) - p_i(z)) / (X-z)
    pub(crate) commitment_to_witness: Commitment<E>,
    /// These are the results of the evaluating each polynomial p_i at the
    /// point `z`.
    pub(crate) evaluated_points: Vec<E::Fr>,
    /// These are the commitments to the p_i polynomials.
    pub(crate) commitments_to_polynomials: Vec<Commitment<E>>,
    _marker: PhantomData<P>,
}

impl<E: PairingEngine<Fr = P::BaseField>, P: TEModelParameters>
    PCAggregateProof<E, P>
{
    /// Initialises an `AggregatedProof` with the commitment to the witness.
    pub(crate) fn with_witness(
        witness: Commitment<E>,
    ) -> PCAggregateProof<E, P> {
        PCAggregateProof {
            commitment_to_witness: witness,
            evaluated_points: Vec::new(),
            commitments_to_polynomials: Vec::new(),
            _marker: PhantomData,
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
    pub(crate) fn flatten(&self, transcript: &mut Transcript) -> KZGProof<E> {
        let challenge =
            transcript.challenge_scalar::<E::Fr>(b"aggregate_witness");
        let powers: Vec<E::Fr> = crate::util::powers_of(
            &challenge,
            self.commitments_to_polynomials.len() - 1,
        );

        let flattened_poly_commitments_iter =
            self.commitments_to_polynomials.iter().zip(powers.iter());
        let flattened_poly_evaluations_iter =
            self.evaluated_points.iter().zip(powers.iter());

        // Flattened polynomial commitments using challenge
        let flattened_poly_commitments = Commitment(
            flattened_poly_commitments_iter
                .map(|(poly_commitment, challenge_power)| {
                    (poly_commitment.0).mul(*challenge_power)
                })
                .sum(),
        );
        //);
        // Flattened evaluation points
        let flattened_poly_evaluations: E::Fr = flattened_poly_evaluations_iter
            .map(|(eval, challenge_power)| *challenge_power * eval)
            .sum();

        KZGProof::<E> {
            random_v: Some(flattened_poly_evaluations),
            w: flattened_poly_commitments.0,
        }
    }
}
