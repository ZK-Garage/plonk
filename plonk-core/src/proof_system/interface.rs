//! Interface to compile, prove, and verify PLONK circuits.

use std::marker::PhantomData;
use ark_poly_commit::{PolynomialCommitment};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, SerializationError, Read, Write};
use merlin::Transcript;

use ark_ec::{ModelParameters, TEModelParameters};
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain,
    UVPolynomial
};
use ark_ff::{self, PrimeField, Field};

use crate::{
    label_polynomial,
    constraint_system::{StandardComposer, Variable}, 
    commitment::HomomorphicCommitment, 
    proof_system::{
        arithmetic,
        ProverKey,
        VerifierKey,
        quotient_poly,
        linearisation_poly,
        proof::Proof,
        Verifier,
    },
    error::{Error, to_pc_error},
    transcript::TranscriptProtocol,
};

fn to_scalars<F, P>(
    circuit: &StandardComposer<F, P>,
    vars: &[Variable]
) -> Vec<F> 
where
    P: ModelParameters<BaseField = F>,
    F: PrimeField,
{
    vars.iter().map(|var| circuit.variables[var]).collect()
}

/// Split `t(X)` poly into 4 n-sized polynomials.
#[allow(clippy::type_complexity)] // NOTE: This is an ok type for internal use.
fn split_tx_poly<F>(
    n: usize,
    t_x: &DensePolynomial<F>,
) -> (
    DensePolynomial<F>,
    DensePolynomial<F>,
    DensePolynomial<F>,
    DensePolynomial<F>,
)
where F: PrimeField
{
    (
        DensePolynomial::from_coefficients_vec(t_x[0..n].to_vec()),
        DensePolynomial::from_coefficients_vec(t_x[n..2 * n].to_vec()),
        DensePolynomial::from_coefficients_vec(t_x[2 * n..3 * n].to_vec()),
        DensePolynomial::from_coefficients_vec(t_x[3 * n..].to_vec()),
    )
}

/// Build PI vector for Proof verifications.
fn build_pi<'a, F>(
    pub_input_values: impl IntoIterator<Item = &'a F>,
    pub_input_pos: &[usize],
    trim_size: usize,
) -> Vec<F>
where
    F: Field,
{
    let mut pi = vec![F::zero(); trim_size];
    pub_input_values
        .into_iter()
        .zip(pub_input_pos.iter().copied())
        .for_each(|(value, pos)| {
            pi[pos] = -*value;
        });
    pi
}

/// Compile a circuit
pub fn compile<F, P, PC>(
    public_parameters: &PC::UniversalParams,
    circuit: &mut StandardComposer<F, P>,
    circuit_size: usize,
) -> Result<(ProvingKey<F, PC>, VerifyingKey<F, PC>), Error>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
    PC: HomomorphicCommitment<F>,
{
    // Setup Public Parameters
    let (committer_key, vk) = PC::trim(
        public_parameters, 
        circuit_size, 
        0, 
        None
    ).map_err(to_pc_error::<F, PC>)?;
    
    // Generate & save `ProverKey` with some random values.
    let mut preprocessed_prover_transcript = Transcript::new(b"CircuitCompilation");
    let prover_key = circuit.preprocess_prover(
        &committer_key, 
        &mut preprocessed_prover_transcript, 
        PhantomData::<PC>,
    )?;

    // Generate & save `VerifierKey` with some random values.
    let mut preprocessed_verifier_transcript = Transcript::new(b"CircuitCompilation");
    let verifier_key = circuit.preprocess_verifier(
        &committer_key, 
        &mut preprocessed_verifier_transcript, 
        PhantomData::<PC>,
    )?;

    let public_input_pos = circuit.pi_positions();

    let verifying_key = VerifyingKey {
        verifier_key: verifier_key,
        public_input_pos: public_input_pos,
        vk: vk,
    };

    let proving_key = ProvingKey {
        prover_key: prover_key,
        committer_key: committer_key,
    };

    Ok((proving_key, verifying_key))
}

/// Struct for proving key
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(bound = ""),
)]
pub struct ProvingKey<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    prover_key: ProverKey<F>,
    committer_key: <PC as PolynomialCommitment<F, DensePolynomial<F>>>::CommitterKey,
}

/// Struct for verifying key
#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(
    Clone(bound = ""),
    Debug(
        bound = "arithmetic::VerifierKey<F,PC>: std::fmt::Debug, PC::Commitment: std::fmt::Debug"
    ),
)]
pub struct VerifyingKey<F, PC>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    verifier_key: VerifierKey<F, PC>,
    public_input_pos: Vec<usize>,
    vk: <PC as PolynomialCommitment<F, DensePolynomial<F>>>::VerifierKey,
}

/// Generate proof
pub fn prove<F, P, PC>(
    proving_key: &ProvingKey<F, PC>,
    circuit: StandardComposer<F, P>,
) -> Result<Proof<F, PC>, Error>
where
    P: TEModelParameters<BaseField = F>,
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
{
    let prover_key = proving_key.prover_key.clone();
    let mut transcript = Transcript::new(b"Test");

    let domain: GeneralEvaluationDomain<F> =
    GeneralEvaluationDomain::new(circuit.circuit_size()).ok_or(Error::InvalidEvalDomainSize {
        log_size_of_group: circuit.circuit_size().trailing_zeros(),
        adicity: <<F as ark_ff::FftField>::FftParams as ark_ff::FftParameters>::TWO_ADICITY,
    })?;
    let n = domain.size();

    // 1. Compute witness Polynomials
    //
    // Convert Variables to scalars padding them to the
    // correct domain size.
    let pad = vec![F::zero(); n - circuit.w_l.len()];
    let w_l_scalar = &[&to_scalars(&circuit, &circuit.w_l)[..], &pad].concat();
    let w_r_scalar = &[&to_scalars(&circuit, &circuit.w_r)[..], &pad].concat();
    let w_o_scalar = &[&to_scalars(&circuit, &circuit.w_o)[..], &pad].concat();
    let w_4_scalar = &[&to_scalars(&circuit, &circuit.w_4)[..], &pad].concat();

    // Witnesses are now in evaluation form, convert them to coefficients
    // so that we may commit to them.
    let w_l_poly =
        DensePolynomial::from_coefficients_vec(domain.ifft(w_l_scalar));
    let w_r_poly =
        DensePolynomial::from_coefficients_vec(domain.ifft(w_r_scalar));
    let w_o_poly =
        DensePolynomial::from_coefficients_vec(domain.ifft(w_o_scalar));
    let w_4_poly =
        DensePolynomial::from_coefficients_vec(domain.ifft(w_4_scalar));

    // Add blinders
    let w_polys = [
        label_polynomial!(w_l_poly),
        label_polynomial!(w_r_poly),
        label_polynomial!(w_o_poly),
        label_polynomial!(w_4_poly),
    ];

    // Commit to witness polynomials.
    let (w_commits, w_rands) = PC::commit(&proving_key.committer_key, w_polys.iter(), None)
        .map_err(to_pc_error::<F, PC>)?;

    // Add witness polynomial commitments to transcript.
    //transcript.append_commitments(&*w_commits, PhantomData::<PC>);
    transcript.append(b"w_l", w_commits[0].commitment());
    transcript.append(b"w_r", w_commits[1].commitment());
    transcript.append(b"w_o", w_commits[2].commitment());
    transcript.append(b"w_4", w_commits[3].commitment());

    // 2. Compute permutation polynomial
    //
    // Compute permutation challenges; `beta` and `gamma`.
    let beta = transcript.challenge_scalar(b"beta");
    transcript.append(b"beta", &beta);
    let gamma = transcript.challenge_scalar(b"gamma");
    transcript.append(b"gamma", &gamma);
    assert!(beta != gamma, "challenges must be different");

    let z_poly = circuit.perm.compute_permutation_poly(
        &domain,
        (w_l_scalar, w_r_scalar, w_o_scalar, w_4_scalar),
        beta,
        gamma,
        (
            &prover_key.permutation.left_sigma.0,
            &prover_key.permutation.right_sigma.0,
            &prover_key.permutation.out_sigma.0,
            &prover_key.permutation.fourth_sigma.0,
        ),
    );

    // Commit to permutation polynomial.
    let (z_poly_commit, _) =
        PC::commit(&proving_key.committer_key, &[label_polynomial!(z_poly)], None)
            .map_err(to_pc_error::<F, PC>)?;

    // Add permutation polynomial commitment to transcript.
    transcript.append(b"z", z_poly_commit[0].commitment());

    // 3. Compute public inputs polynomial.
    let pi_poly = DensePolynomial::from_coefficients_vec(
        domain.ifft(&circuit.construct_dense_pi_vec()),
    );

    // 4. Compute quotient polynomial
    //
    // Compute quotient challenge; `alpha`, and gate-specific separation
    // challenges.
    let alpha = transcript.challenge_scalar(b"alpha");
    let range_sep_challenge =
        transcript.challenge_scalar(b"range separation challenge");
    let logic_sep_challenge =
        transcript.challenge_scalar(b"logic separation challenge");
    let fixed_base_sep_challenge =
        transcript.challenge_scalar(b"fixed base separation challenge");
    let var_base_sep_challenge =
        transcript.challenge_scalar(b"variable base separation challenge");

    let t_poly = quotient_poly::compute::<F, P>(
        &domain,
        &prover_key,
        &z_poly,
        &w_l_poly,
        &w_r_poly,
        &w_o_poly,
        &w_4_poly,
        &pi_poly,
        &alpha,
        &beta,
        &gamma,
        &range_sep_challenge,
        &logic_sep_challenge,
        &fixed_base_sep_challenge,
        &var_base_sep_challenge,
    )?;

    let (t_1_poly, t_2_poly, t_3_poly, t_4_poly) =
        split_tx_poly(n, &t_poly);

    // Commit to splitted quotient polynomial
    let (t_commits, _) = PC::commit(
        &proving_key.committer_key,
        &[
            label_polynomial!(t_1_poly),
            label_polynomial!(t_2_poly),
            label_polynomial!(t_3_poly),
            label_polynomial!(t_4_poly),
        ],
        None,
    )
    .map_err(to_pc_error::<F, PC>)?;

    // Add quotient polynomial commitments to transcript
    transcript.append(b"t_1", t_commits[0].commitment());
    transcript.append(b"t_2", t_commits[1].commitment());
    transcript.append(b"t_3", t_commits[2].commitment());
    transcript.append(b"t_4", t_commits[3].commitment());

    // 4. Compute linearisation polynomial
    //
    // Compute evaluation challenge; `z`.
    let z_challenge = transcript.challenge_scalar(b"z");

    let (lin_poly, evaluations) = linearisation_poly::compute::<F, P>(
        &domain,
        &prover_key,
        &alpha,
        &beta,
        &gamma,
        &range_sep_challenge,
        &logic_sep_challenge,
        &fixed_base_sep_challenge,
        &var_base_sep_challenge,
        &z_challenge,
        &w_l_poly,
        &w_r_poly,
        &w_o_poly,
        &w_4_poly,
        &t_1_poly,
        &t_2_poly,
        &t_3_poly,
        &t_4_poly,
        &z_poly,
    )?;

    // Add evaluations to transcript.
    // First wire evals
    transcript.append(b"a_eval", &evaluations.wire_evals.a_eval);
    transcript.append(b"b_eval", &evaluations.wire_evals.b_eval);
    transcript.append(b"c_eval", &evaluations.wire_evals.c_eval);
    transcript.append(b"d_eval", &evaluations.wire_evals.d_eval);

    // Second permutation evals
    transcript
        .append(b"left_sig_eval", &evaluations.perm_evals.left_sigma_eval);
    transcript.append(
        b"right_sig_eval",
        &evaluations.perm_evals.right_sigma_eval,
    );
    transcript
        .append(b"out_sig_eval", &evaluations.perm_evals.out_sigma_eval);
    transcript
        .append(b"perm_eval", &evaluations.perm_evals.permutation_eval);

    // Third, all evals needed for custom gates
    evaluations
        .custom_evals
        .vals
        .iter()
        .for_each(|(label, eval)| {
            let static_label = Box::leak(label.to_owned().into_boxed_str());
            transcript.append(static_label.as_bytes(), eval);
        });

    // 5. Compute Openings using KZG10
    //
    // We merge the quotient polynomial using the `z_challenge` so the SRS
    // is linear in the circuit size `n`

    // Compute aggregate witness to polynomials evaluated at the evaluation
    // challenge `z`
    let aw_challenge: F = transcript.challenge_scalar(b"aggregate_witness");

    let aw_polys = [
        label_polynomial!(lin_poly),
        label_polynomial!(prover_key.permutation.left_sigma.0.clone()),
        label_polynomial!(prover_key.permutation.right_sigma.0.clone()),
        label_polynomial!(prover_key.permutation.out_sigma.0.clone()),
    ];

    let (aw_commits, aw_rands) = PC::commit(&proving_key.committer_key, &aw_polys, None)
        .map_err(to_pc_error::<F, PC>)?;

    let aw_opening = PC::open(
        &proving_key.committer_key,
        aw_polys.iter().chain(w_polys.iter()),
        aw_commits.iter().chain(w_commits.iter()),
        &z_challenge,
        aw_challenge,
        aw_rands.iter().chain(w_rands.iter()),
        None,
    )
    .map_err(to_pc_error::<F, PC>)?;

    let saw_challenge: F =
        transcript.challenge_scalar(b"aggregate_witness");

    let saw_polys = [
        label_polynomial!(z_poly),
        label_polynomial!(w_l_poly),
        label_polynomial!(w_r_poly),
        label_polynomial!(w_4_poly),
    ];

    let (saw_commits, saw_rands) = PC::commit(&proving_key.committer_key, &saw_polys, None)
        .map_err(to_pc_error::<F, PC>)?;

    let saw_opening = PC::open(
        &proving_key.committer_key,
        &saw_polys,
        &saw_commits,
        &(z_challenge * domain.element(1)),
        saw_challenge,
        &saw_rands,
        None,
    )
    .map_err(to_pc_error::<F, PC>)?;

    Ok(Proof {
        a_comm: w_commits[0].commitment().clone(),
        b_comm: w_commits[1].commitment().clone(),
        c_comm: w_commits[2].commitment().clone(),
        d_comm: w_commits[3].commitment().clone(),
        z_comm: saw_commits[0].commitment().clone(),
        t_1_comm: t_commits[0].commitment().clone(),
        t_2_comm: t_commits[1].commitment().clone(),
        t_3_comm: t_commits[2].commitment().clone(),
        t_4_comm: t_commits[3].commitment().clone(),
        aw_opening,
        saw_opening,
        evaluations,
    })
}

/// Verify a proof
pub fn verify<F, P, PC>(
    verifying_key: &VerifyingKey<F, PC>,
    public_input: &[F],
    proof: &Proof<F, PC>
) -> Result<bool, Error>
where
    F: PrimeField,
    PC: HomomorphicCommitment<F>,
    P: TEModelParameters<BaseField = F>,
{
    let mut verifier: Verifier<F, P, PC> = Verifier::new(b"Test");
    let padded_circuit_size = verifying_key.verifier_key.padded_circuit_size();
    verifier.verifier_key = Some(verifying_key.verifier_key.clone());

    let res = verifier.verify(
        proof,
        &verifying_key.vk,
        build_pi(public_input, &verifying_key.public_input_pos, padded_circuit_size)
            .as_slice(),
    );

    match res {
        Ok(()) => Ok(true),
        Err(Error::ProofVerificationError) => Ok(false),
        Err(e) => panic!("{:?}", e),
    }
}
