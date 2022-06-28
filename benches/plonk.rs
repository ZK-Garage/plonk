// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! PLONK Benchmarks

use ark_poly_commit::PolynomialCommitment;
use core::marker::PhantomData;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use plonk::{parameters::test::*, prelude::*};
use rand_core::OsRng;

/// Benchmark Circuit
#[derive(derivative::Derivative)]
#[derivative(Debug, Default)]
pub struct BenchCircuit<P> {
    /// Circuit Size
    size: usize,

    /// Field and parameters
    _phantom: PhantomData<P>,
}

impl<P> BenchCircuit<P> {
    /// Builds a new circuit with a constraint count of `2^degree`.
    #[inline]
    pub fn new(degree: usize) -> Self {
        Self {
            size: 1 << degree,
            _phantom: PhantomData::<P>,
        }
    }
}

impl<P> Circuit<P> for BenchCircuit<P>
where
    P: CircuitParameters,
{
    const CIRCUIT_ID: [u8; 32] = [0xff; 32];

    #[inline]
    fn gadget(
        &mut self,
        composer: &mut StandardComposer<P>,
    ) -> Result<(), Error> {
        composer.add_dummy_lookup_table();
        while composer.circuit_bound() < self.size - 1 {
            composer.add_dummy_constraints();
        }
        Ok(())
    }

    #[inline]
    fn padded_circuit_size(&self) -> usize {
        self.size
    }
}

fn kzg10_benchmarks(c: &mut Criterion) {
    constraint_system_benchmark::<Bls12_381_KZG>("KZG10", c);
}

fn ipa_benchmarks(c: &mut Criterion) {
    constraint_system_benchmark::<Bls12_377_IPA>("IPA", c);
}

/// Generates full benchmark suite for compiling, proving, and verifying.
fn constraint_system_benchmark<P>(name: &str, c: &mut Criterion)
where
    P: CircuitParameters,
{
    let label = b"ark".as_slice();

    const MINIMUM_DEGREE: usize = 5;
    const MAXIMUM_DEGREE: usize = 19;

    let pp =
        P::PolynomialCommitment::setup(1 << MAXIMUM_DEGREE, None, &mut OsRng)
            .expect("Unable to sample public parameters.");

    let mut compiling_benchmarks =
        c.benchmark_group(format!("{0}/compile", name));
    for degree in MINIMUM_DEGREE..MAXIMUM_DEGREE {
        let mut circuit = BenchCircuit::<P>::new(degree);
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
        let mut circuit = BenchCircuit::<P>::new(degree);
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
        let mut circuit = BenchCircuit::<P>::new(degree);
        let (pk_p, vk) =
            circuit.compile(&pp).expect("Unable to compile circuit.");
        let (proof, pi) = circuit.gen_proof(&pp, pk_p.clone(), &label).unwrap();
        verifying_benchmarks.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    plonk::circuit::verify_proof::<P>(
                        &pp,
                        vk.clone(),
                        &proof,
                        &pi,
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
    targets = kzg10_benchmarks, ipa_benchmarks
}
criterion_main!(plonk);
