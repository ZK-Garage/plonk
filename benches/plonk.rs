use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
use ark_ec::{PairingEngine, TEModelParameters};
use ark_ed_on_bls12_381::EdwardsParameters;
use ark_plonk::prelude::*;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::kzg10::KZG10;
use criterion::{
    black_box, criterion_group, criterion_main, BenchmarkId, Criterion,
};
use rand_core::OsRng;
use std::marker::PhantomData;

#[derive(Debug)]
pub struct BenchCircuit<
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
> {
    degree: usize,
    _marker: PhantomData<P>,
    _marker2: PhantomData<E>,
}

impl<E: PairingEngine, P: TEModelParameters<BaseField = E::Fr>> Default
    for BenchCircuit<E, P>
{
    fn default() -> Self {
        Self {
            degree: 1 << 2,
            _marker: PhantomData,
            _marker2: PhantomData,
        }
    }
}

impl<
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
        T: Into<usize>,
    > From<T> for BenchCircuit<E, P>
{
    fn from(degree: T) -> Self {
        Self {
            degree: 1 << degree.into(),
            _marker: PhantomData,
            _marker2: PhantomData,
        }
    }
}

impl<E: PairingEngine, P: TEModelParameters<BaseField = E::Fr>> Circuit<E, P>
    for BenchCircuit<E, P>
{
    const CIRCUIT_ID: [u8; 32] = [0xff; 32];

    fn gadget(
        &mut self,
        composer: &mut StandardComposer<E, P>,
    ) -> Result<(), Error> {
        while composer.circuit_size() < self.degree - 1 {
            composer.add_dummy_constraints();
        }

        Ok(())
    }

    fn padded_circuit_size(&self) -> usize {
        self.degree
    }
}

impl<E: PairingEngine, P: TEModelParameters<BaseField = E::Fr>>
    BenchCircuit<E, P>
{
    pub fn degree(&self) -> usize {
        self.degree
    }
}

fn constraint_system_benchmark(c: &mut Criterion) {
    let initial_degree = 5usize;
    let final_degree = 18usize;
    let label = b"ark".as_slice();

    // Generate CRS
    let pp = KZG10::<Bls12_381, DensePolynomial<BlsScalar>>::setup(
        1 << 19,
        false,
        &mut OsRng,
    )
    .expect("UP generation err");

    for degree in initial_degree..=final_degree {
        // Gen circuit instance with this round degree
        let mut circuit =
            BenchCircuit::<Bls12_381, EdwardsParameters>::from(degree);
        // Compile the circuit
        let (pk_p, verifier_data) =
            circuit.compile(&pp).expect("Compilation err");
        let description = format!("Prove 2^{} = {} gates", degree, 1 << degree);

        c.bench_with_input(
            BenchmarkId::new(description.as_str(), degree),
            &degree,
            |b, _degree| {
                b.iter(|| {
                    circuit
                        .gen_proof(
                            black_box(&pp),
                            black_box(pk_p.clone()),
                            black_box(&label),
                        )
                        .unwrap()
                })
            },
        );

        let description =
            format!("Verify 2^{} = {} gates", degree, 1 << degree);
        let proof = circuit
            .gen_proof(
                black_box(&pp),
                black_box(pk_p.clone()),
                black_box(&label),
            )
            .unwrap();

        c.bench_with_input(
            BenchmarkId::new(description.as_str(), degree),
            &degree,
            |b, _degree| {
                b.iter(|| {
                    ark_plonk::circuit::verify_proof(
                        black_box(&pp),
                        black_box(verifier_data.clone().key()),
                        black_box(&proof),
                        &[],
                        black_box(&verifier_data.clone().pi_pos()),
                        black_box(&label),
                    )
                    .expect("Failed to verify bench circuit!");
                })
            },
        );
    }
}

criterion_group! {
    name = ark_plonk;
    config = Criterion::default().sample_size(10);
    targets = constraint_system_benchmark
}
criterion_main!(ark_plonk);
