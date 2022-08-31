// Licensed under the Apache License, Version 2.0 <LICENSE-APACHE
// or https://www.apache.org/licenses/LICENSE-2.0> or the MIT license
// <LICENSE-MIT or https://opensource.org/licenses/MIT>, at your
// option. This file may not be copied, modified, or distributed
// except according to those terms.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! PLONK Benchmarks

use ark_bls12_377::Bls12_377;
use ark_bls12_381::Bls12_381;
use ark_ec::{PairingEngine, TEModelParameters};
use ark_ed_on_bls12_381::EdwardsParameters;
use ark_ff::{FftField, PrimeField};
use blake2;
use core::marker::PhantomData;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use plonk::commitment::{HomomorphicCommitment, IPA, KZG10};
use plonk::prelude::*;
use rand_core::OsRng;

/// Benchmark Circuit
#[derive(derivative::Derivative)]
#[derivative(Debug, Default)]
pub struct BenchCircuit<F, P> {
    /// Circuit Size
    size: usize,

    /// Field and parameters
    _phantom: PhantomData<(F, P)>,
}

impl<F, P> BenchCircuit<F, P> {
    /// Builds a new circuit with a constraint count of `2^degree`.
    #[inline]
    pub fn new(degree: usize) -> Self {
        Self {
            size: 1 << degree,
            _phantom: PhantomData::<(F, P)>,
        }
    }
}

impl<F, P> Circuit<F, P> for BenchCircuit<F, P>
where
    F: FftField + PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    const CIRCUIT_ID: [u8; 32] = [0xff; 32];

    #[inline]
    fn gadget(
        &mut self,
        composer: &mut StandardComposer<F, P>,
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
    constraint_system_benchmark::<
        <Bls12_381 as PairingEngine>::Fr,
        EdwardsParameters,
        KZG10<Bls12_381>,
    >("KZG10", c);
}

fn ipa_benchmarks(c: &mut Criterion) {
    constraint_system_benchmark::<
        <Bls12_377 as PairingEngine>::Fr,
        ark_ed_on_bls12_377::EdwardsParameters,
        IPA<<Bls12_377 as PairingEngine>::G1Affine, blake2::Blake2b>,
    >("IPA", c);
}

/// Generates full benchmark suite for compiling, proving, and verifying.
fn constraint_system_benchmark<F, P, HC>(name: &str, c: &mut Criterion)
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
    HC: HomomorphicCommitment<F>,
{
    let label = b"ark".as_slice();

    const MINIMUM_DEGREE: usize = 5;
    const MAXIMUM_DEGREE: usize = 19;

    let pp = HC::setup(1 << MAXIMUM_DEGREE, None, &mut OsRng)
        .expect("Unable to sample public parameters.");

    let mut compiling_benchmarks =
        c.benchmark_group(format!("{0}/compile", name));
    for degree in MINIMUM_DEGREE..MAXIMUM_DEGREE {
        let mut circuit = BenchCircuit::<F, P>::new(degree);
        compiling_benchmarks.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    circuit
                        .compile::<HC>(&pp)
                        .expect("Unable to compile circuit.")
                })
            },
        );
    }
    compiling_benchmarks.finish();

    let mut proving_benchmarks = c.benchmark_group("prove");
    for degree in MINIMUM_DEGREE..MAXIMUM_DEGREE {
        let mut circuit = BenchCircuit::<F, P>::new(degree);
        let (pk_p, _) = circuit
            .compile::<HC>(&pp)
            .expect("Unable to compile circuit.");
        proving_benchmarks.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    circuit.gen_proof::<HC>(&pp, pk_p.clone(), &label).unwrap()
                })
            },
        );
    }
    proving_benchmarks.finish();

    let mut verifying_benchmarks = c.benchmark_group("verify");
    for degree in MINIMUM_DEGREE..MAXIMUM_DEGREE {
        let mut circuit = BenchCircuit::<F, P>::new(degree);
        let (pk_p, (vk, _pi_pos)) =
            circuit.compile(&pp).expect("Unable to compile circuit.");
        let (proof, pi) =
            circuit.gen_proof::<HC>(&pp, pk_p.clone(), &label).unwrap();
        verifying_benchmarks.bench_with_input(
            BenchmarkId::from_parameter(degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    plonk::circuit::verify_proof::<F, P, HC>(
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
