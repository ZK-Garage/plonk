# PLONK
![Build Status](https://github.com/rust-zkp/ark-plonk/workflows/Continuous%20integration/badge.svg)
[![Repository](https://img.shields.io/badge/github-plonk-blueviolet?logo=github)](https://github.com/rust-zkp/ark-plonk)
[![Documentation](https://img.shields.io/badge/docs-plonk-blue?logo=rust)](https://docs.rs/plonk/)


_This is a pure Rust implementation of the PLONK zk proving system_


## About
Initial implementation created by [Kev](https://github.com/kevaundray), [Carlos](https://github.com/CPerezz) and [Luke](https://github.com/LukePearson1) at Dusk Network.
Redesigned by the [rust zkp](https://github.com/rust-zkp) team to have a backend which is compatible with the [arkworks](https://github.com/arkworks-rs) suite. This allows us to leverage the multitude of curves
and optimised algebra present in various arkworks repositories.

## Usage
```rust
use core::marker::PhantomData;
use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
use ark_ec::twisted_edwards_extended::GroupAffine;
use ark_ec::{AffineCurve, PairingEngine, TEModelParameters};
use ark_ed_on_bls12_381::{
    EdwardsAffine as JubjubAffine, EdwardsParameters as JubjubParameters,
    EdwardsProjective as JubjubProjective, Fr as JubjubScalar,
};
use ark_ff::{BigInteger, PrimeField};
use ark_plonk::circuit::{self, Circuit, PublicInputValue};
use ark_plonk::prelude::*;
use num_traits::{One, Zero};
use rand_core::OsRng;

// Implement a circuit that checks:
// 1) a + b = c where C is a PI
// 2) a <= 2^6
// 3) b <= 2^5
// 4) a * b = d where D is a PI
// 5) JubJub::GENERATOR * e(JubJubScalar) = f where F is a Public Input
#[derive(Debug, Default)]
pub struct TestCircuit<
    E: PairingEngine,
    P: TEModelParameters<BaseField = E::Fr>,
> {
    a: E::Fr,
    b: E::Fr,
    c: E::Fr,
    d: E::Fr,
    e: P::ScalarField,
    f: GroupAffine<P>,
}
impl<
        E: PairingEngine,
        P: TEModelParameters<BaseField = E::Fr>,
    > Circuit<E, P> for TestCircuit<E, P>
{
    const CIRCUIT_ID: [u8; 32] = [0xff; 32];
    fn gadget(
        &mut self,
        composer: &mut StandardComposer<E, P>,
    ) -> Result<(), Error> {
        let a = composer.add_input(self.a);
        let b = composer.add_input(self.b);
        // Make first constraint a + b = c
        let add_result = composer.add(
          (E::Fr::one(), a),
          (E::Fr::one(), b),
          E::Fr::zero(),
          Some(-self.c),
        );
	composer.assert_equal(add_result, composer.zero_var());

        // Check that a and b are in range
        composer.range_gate(a, 1 << 6);
        composer.range_gate(b, 1 << 5);
        // Make second constraint a * b = d
        let mul_result = composer.mul(E::Fr::one(), a, b, E::Fr::zero(), Some(-self.d));
        composer.assert_equal(mul_result, composer.zero_var());

        let e_repr = self.e.into_repr().to_bytes_le();
        let e = composer.add_input(E::Fr::from_le_bytes_mod_order(&e_repr));
        let (x, y) = P::AFFINE_GENERATOR_COEFFS;
        let generator = GroupAffine::new(x, y);
        let scalar_mul_result = composer.fixed_base_scalar_mul(e, generator);
        // Apply the constrain
        composer.assert_equal_public_point(scalar_mul_result, self.f);
        Ok(())
    }
    fn padded_circuit_size(&self) -> usize {
        1 << 11
    }
}

// Now let's use the Circuit we've just implemented!
fn main()-> Result<(), Error> {
    let pp: PublicParameters<Bls12_381> = KZG10::<Bls12_381,DensePolynomial<BlsScalar>,>::setup(
          1 << 12, false, &mut OsRng
    )?;
    // Initialize the circuit
    let mut circuit: TestCircuit::<
        Bls12_381,
        JubjubParameters,
    > = TestCircuit::default();
    // Compile the circuit
    let (pk, vd) = circuit.compile(&pp).unwrap();
    // Generator
    let (x, y) = JubJubParameters::AFFINE_GENERATOR_COEFFS;
    let generator = JubJubAffine::new(x, y);
    let point_f_pi: JubJubAffine = AffineCurve::mul(
      &generator,
      JubJubScalar::from(2u64).into_repr(),
    ).into_affine();
    let proof = {
        let mut circuit = TestCircuit {
            a: BlsScalar::from(20u64),
            b: BlsScalar::from(5u64),
            c: BlsScalar::from(25u64),
            d: BlsScalar::from(100u64),
            e: JubJubScalar::from(2u64),
            f: point_f_pi,
        };
        circuit.gen_proof(&pp, pk, b"Test").unwrap()
    };

    let public_inputs: Vec<PublicInputValue<BlsScalar, JubjubParameters>> = vec![
        BlsScalar::from(25u64).into_pi(),
        BlsScalar::from(100u64).into_pi(),
        point_f_pi.into_pi()
    ];

    circuit::verify_proof(
        &pp,
        *vd.key(),
        &proof,
        &public_inputs,
        &vd.pi_pos(),
        b"Test",
    )
    .unwrap();
}
```

### Features

This crate includes a variety of features which will briefly be explained below:
- `parallel`: Enables `rayon` and other parallelisation primitives to be used and speed up some of the algorithms used
by the crate and it's dependencies.
- `asm`: Enables inline-assembly implementations for some of the internal algorithms and primitives used by the `arkworks` dependencies of the crate.
- `trace`: Enables the Circuit debugger tooling. This is essentially the capability of using the
`StandardComposer::check_circuit_satisfied` function. The function will output information about each circuit gate until one of the gates does not satisfy the equation, or there are no more gates. If there is an unsatisfied gate
equation, the function will panic and return the gate number.
- `trace-print`: Goes a step further than `trace` and prints each `gate` component data, giving a clear overview of all the
values which make up the circuit that we're constructing.
__The recommended method is to derive the std output, and the std error, and then place them in text file
  which can be used to efficiently analyse the gates.__



## Documentation

There are two main types of documentation in this repository:

- **Crate documentation**. This provides info about all of the functions that the library provides, as well
  as the documentation regarding the data structures that it exports. To check this, please feel free to go to
  the [documentation page](https://docs.rs/ark-plonk/) or run `make doc` or `make doc-internal`.

- **Notes**. This is a specific subset of documentation which explains the key mathematical concepts
  of PLONK and how they work with mathematical demonstrations. To check it, run `make doc` and open the resulting docs,
  which will be located under `/target` with your browser.

## Performance

Benches taken running: `RUSTFLAGS='-C target-cpu=native' cargo bench`
with an `Intel(R) Core(TM) i9-10885H`
```
Prove 2^5 = 32 gates/5  time:   [9.4230 ms 9.5398 ms 9.6632 ms]                                  
Verify 2^5 = 32 gates/5 time:   [4.1958 ms 4.2881 ms 4.4535 ms]                                    

Prove 2^6 = 64 gates/6  time:   [12.896 ms 13.013 ms 13.141 ms]                                  
Verify 2^6 = 64 gates/6 time:   [4.2475 ms 4.2781 ms 4.3019 ms]                                    

Prove 2^7 = 128 gates/7 time:   [17.836 ms 18.137 ms 18.349 ms]                                   
Verify 2^7 = 128 gates/7 time:   [4.2062 ms 4.2973 ms 4.3528 ms]

Prove 2^8 = 256 gates/8 time:   [29.553 ms 29.914 ms 30.252 ms]                                   
Verify 2^8 = 256 gates/8 time:   [4.1951 ms 4.2593 ms 4.3262 ms]

Prove 2^9 = 512 gates/9 time:   [49.285 ms 50.221 ms 51.109 ms]                                   
Verify 2^9 = 512 gates/9 time:   [4.2094 ms 4.3023 ms 4.3846 ms]

Prove 2^10 = 1024 gates/10 time:   [68.426 ms 68.704 ms 68.854 ms]
Verify 2^10 = 1024 gates/10 time:   [4.1785 ms 4.2228 ms 4.2789 ms]

Prove 2^11 = 2048 gates/11 time:   [127.20 ms 127.49 ms 127.70 ms]
Verify 2^11 = 2048 gates/11 time:   [4.1177 ms 4.1379 ms 4.1461 ms]

Prove 2^12 = 4096 gates/12 time:   [245.00 ms 245.48 ms 245.92 ms]
Verify 2^12 = 4096 gates/12 time:   [4.1342 ms 4.1467 ms 4.1632 ms]

Prove 2^13 = 8192 gates/13 time:   [438.23 ms 440.64 ms 442.87 ms]
Verify 2^13 = 8192 gates/13 time:   [4.1329 ms 4.1770 ms 4.2117 ms]

Prove 2^14 = 16384 gates/14 time:   [865.00 ms 869.66 ms 873.92 ms]
Verify 2^14 = 16384 gates/14 time:   [4.1817 ms 4.1870 ms 4.1921 ms]

Prove 2^15 = 32768 gates/15 time:   [1.7644 s 1.7712 s 1.7778 s]
Verify 2^15 = 32768 gates/15 time:   [4.3285 ms 4.3390 ms 4.3591 ms]

Prove 2^16 = 65536 gates/16 time:   [3.4207 s 3.4499 s 3.4776 s]
Verify 2^16 = 65536 gates/16 time:   [4.4808 ms 4.5020 ms 4.5267 ms]

Prove 2^17 = 131072 gates/17 time:   [6.6989 s 6.7577 s 6.8191 s]
Verify 2^17 = 131072 gates/17 time:   [5.1339 ms 5.1572 ms 5.1806 ms]

Prove 2^18 = 262144 gates/18 time:   [13.603 s 13.704 s 13.816 s]
Verify 2^18 = 262144 gates/18 time:   [6.7577 ms 6.8124 ms 6.8925 ms]
```

## Acknowledgements

- Reference [implementation](https://github.com/AztecProtocol/barretenberg) by Aztec Protocol
- Initial [implementation](https://github.com/kobigurk/plonk/tree/kobigurk/port_to_zexe) of PLONK with arkworks backend was done years before this lib existed by Kobi Gurkan

## Licensing

This code is licensed under Mozilla Public License Version 2.0 (MPL-2.0). Please see [LICENSE](https://github.com/rust-zkp/ark-plonk/blob/master/LICENSE) for further info.

## Contributing

- If you want to contribute to this repository/project please, check [CONTRIBUTING.md](https://github.com/rust-zkp/ark-plonk/blob/master/CONTRIBUTING.md)
- If you want to report a bug or request a new feature addition, please open an issue on this repository.
