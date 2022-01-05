// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! PLONK Benchmarks

use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
use ark_ec::{PairingEngine, TEModelParameters};
use ark_ed_on_bls12_381::EdwardsParameters;
use ark_poly::univariate::DensePolynomial;
use ark_poly_commit::kzg10::KZG10;
use core::marker::PhantomData;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use plonk::prelude::*;
use rand_core::OsRng;

/// Benchmark Circuit
#[derive(derivative::Derivative)]
#[derivative(Debug, Default)]
pub struct BenchCircuit<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    /// Circuit Size
    size: usize,

    /// Type Parameter Marker
    __: PhantomData<(E, P)>,
}

impl<E, P> BenchCircuit<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    /// Builds a new circuit with a constraint count of `2^degree`.
    #[inline]
    pub fn new(degree: usize) -> Self {
        Self {
            size: 1 << degree,
            __: PhantomData,
        }
    }
}

impl<E, P> Circuit<E, P> for BenchCircuit<E, P>
where
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
{
    const CIRCUIT_ID: [u8; 32] = [0xff; 32];

    #[inline]
    fn gadget(
        &mut self,
        composer: &mut StandardComposer<E, P>,
    ) -> Result<(), Error> {
        while composer.circuit_size() < self.size - 1 {
            composer.add_dummy_constraints();
        }
        Ok(())
    }

    #[inline]
    fn padded_circuit_size(&self) -> usize {
        self.size
    }
}

/// Generates full benchmark suite for compiling, proving, and verifying.
fn constraint_system_benchmark(c: &mut Criterion) {
    let label = b"ark".as_slice();

    const MINIMUM_DEGREE: usize = 5;
    const MAXIMUM_DEGREE: usize = 19;

    let pp = KZG10::<Bls12_381, DensePolynomial<BlsScalar>>::setup(
        // +1 per wire, +2 for the permutation poly
        1 << MAXIMUM_DEGREE + 6,
        false,
        &mut OsRng,
    )
    .expect("Unable to sample public parameters.");

    let mut compiling_benchmarks = c.benchmark_group("compile");
    for degree in MINIMUM_DEGREE..MAXIMUM_DEGREE {
        let mut circuit = BenchCircuit::<_, EdwardsParameters>::new(degree);
        compiling_benchmarks.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    circuit.compile(&pp).expect("Unable to compile circuit.")
                })
            },
        );
    }
    compiling_benchmarks.finish();

    let mut proving_benchmarks = c.benchmark_group("prove");
    for degree in MINIMUM_DEGREE..MAXIMUM_DEGREE {
        let mut circuit = BenchCircuit::<_, EdwardsParameters>::new(degree);
        let (pk_p, _) =
            circuit.compile(&pp).expect("Unable to compile circuit.");
        proving_benchmarks.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, _| {
                b.iter(|| circuit.gen_proof(&pp, pk_p.clone(), &label).unwrap())
            },
        );
    }
    proving_benchmarks.finish();

    let mut verifying_benchmarks = c.benchmark_group("verify");
    for degree in MINIMUM_DEGREE..MAXIMUM_DEGREE {
        let mut circuit = BenchCircuit::<_, EdwardsParameters>::new(degree);
        let (pk_p, verifier_data) =
            circuit.compile(&pp).expect("Unable to compile circuit.");
        let proof = circuit.gen_proof(&pp, pk_p.clone(), &label).unwrap();
        let VerifierData { key, pi_pos } = verifier_data;
        verifying_benchmarks.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    plonk::circuit::verify_proof(
                        &pp,
                        key.clone(),
                        &proof,
                        &[],
                        &pi_pos,
                        &label,
                    )
                    .expect("Unable to verify benchmark circuit.");
                })
            },
        );
    }
    verifying_benchmarks.finish();
}

criterion_group! {
    name = plonk;
    config = Criterion::default().sample_size(10);
    targets = constraint_system_benchmark
}
criterion_main!(plonk);
