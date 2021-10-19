# PLONK 
![Build Status](https://github.com/rust-zkp/ark-plonk/workflows/Continuous%20integration/badge.svg)
[![Repository](https://img.shields.io/badge/github-plonk-blueviolet?logo=github)](https://github.com/rust-zkp/ark-plonk)
[![Documentation](https://img.shields.io/badge/docs-plonk-blue?logo=rust)](https://docs.rs/plonk/)


_This is a pure Rust implementation of the PLONK zk proving system_

## Usage

```rust
use ark_plonk::prelude::*;
use ark_ec::bls12::Bls12;
use rand_core::OsRng;


// Implement a circuit that checks:
// 1) a + b = c where C is a PI
// 2) a <= 2^6
// 3) b <= 2^5
// 4) a * b = d where D is a PI
// 5) JubJub::GENERATOR * e(JubJubScalar) = f where F is a Public Input
#[derive(Debug, Default)]
pub struct TestCircuit {
    a: Fr,
    b: Fr,
    c: Fr,
    d: Fr,
    e: JubJubScalar,
    f: JubJubAffine,
}
///XXX: TODO: Swap out JubJub for arkworks embedded curves
impl Circuit for TestCircuit {
    const CIRCUIT_ID: [u8; 32] = [0xff; 32];
    fn gadget(
        &mut self,
        composer: &mut StandardComposer,
    ) -> Result<(), Error> {
        let a = composer.add_input(self.a);
        let b = composer.add_input(self.b);
        // Make first constraint a + b = c
        composer.poly_gate(
            a,
            b,
            composer.zero_var(),
            Bls12::Fr::zero(),
            Bls12::Fr::one(),
            Bls12::Fr::one(),
            Bls12::Fr::zero(),
            Bls12::Fr::zero(),
            Some(-self.c),
        );
        // Check that a and b are in range
        composer.range_gate(a, 1 << 6);
        composer.range_gate(b, 1 << 5);
        // Make second constraint a * b = d
        composer.poly_gate(
            a,
            b,
            composer.zero_var(),
            Bls12::Fr::one(),
            Bls12::Fr::zero(),
            Bls12::Fr::zero(),
            Bls12::Fr::one(),
            Bls12::Fr::zero(),
            Some(-self.d),
        );

        let e = composer.add_input(self.e.into());
        let scalar_mul_result = composer
            .fixed_base_scalar_mul(e, dusk_jubjub::GENERATOR_EXTENDED);
        // Apply the constrain
        composer.assert_equal_public_point(scalar_mul_result, self.f);
        Ok(())
    }
    fn padded_circuit_size(&self) -> usize {
        1 << 11
    }
}

// Now let's use the Circuit we've just implemented!

let pp = PublicParameters::setup(1 << 12, &mut OsRng).unwrap();
// Initialize the circuit
let mut circuit = TestCircuit::default();
// Compile the circuit
let (pk, vd) = circuit.compile(&pp).unwrap();
// Prover POV
let proof = {
    let mut circuit = TestCircuit {
        a: Bls12::Fr::from(20u64),
        b: Bls12::Fr::from(5u64),
        c: Bls12::Fr::from(25u64),
        d: Bls12::Fr::from(100u64),
        e: JubJubScalar::from(2u64),
        f: JubJubAffine::from(
            dusk_jubjub::GENERATOR_EXTENDED * JubJubScalar::from(2u64),
        ),
    };
    circuit.gen_proof(&pp, &pk, b"Test").unwrap()
};
// Verifier POV
let public_inputs: Vec<PublicInputValue> = vec![
    Bls12::Fr::from(25u64).into(),
    Bls12::Fr::from(100u64).into(),
    JubJubAffine::from(
        dusk_jubjub::GENERATOR_EXTENDED * JubJubScalar::from(2u64),
    )
    .into(),
];
circuit::verify_proof(
    &pp,
    &vd.key(),
    &proof,
    &public_inputs,
    &vd.pi_pos(),
    b"Test",
).unwrap();
```

### Features

This crate includes a variety of features which will briefly be explained below:
- `alloc`: Enables the usage of an allocator and with it the capability of performing `Proof` constructions and 
  verifications. Without this feature it **IS NOT** possible to prove or verify anything. 
  Its absence only makes `ark-plonk` export certain fixed-size data structures such as `Proof` which can be useful in no_std envoirments where we don't have allocators either.
- `std`: Enables `std` usage as well as `rayon` parallelisation in some proving and verifying ops. 
- `trace`: Enables the Circuit debugger tooling. This is essentially the capability of using the 
  `StandardComposer::check_circuit_satisfied` function. The function will output information about each circuit gate until 
  one of the gates does not satisfy the equation, or there are no more gates. If there is an unsatisfied gate 
  equation, the function will panic and return the gate number.
- `trace-print`: Goes a step further than `trace` and prints each `gate` component data, giving a clear overview of all the 
  values which make up the circuit that we're constructing. 
  __The recommended method is to derive the std output, and the std error, and then place them in text file 
    which can be used to efficiently analyse the gates.__



## Documentation

There are two main types of documentation in this repository:

- **Crate documentation**. This provides info about all of the functions that the library provides, as well
  as the documentation regarding the data structures that it exports. To check this, please feel free to go to
  the [documentation page](https://docs.rs/dusk-plonk/) or run `make doc` or `make doc-internal`.

- **Notes**. This is a specific subset of documentation which explains the key mathematical concepts
  of PLONK and how they work with mathematical demonstrations. To check it, run `make doc` and open the resulting docs,
  which will be located under `/target` with your browser.

## Performance
TODO

## Acknowledgements

- Reference [implementation](https://github.com/AztecProtocol/barretenberg) by Aztec Protocol
- Initial [implementation](https://github.com/kobigurk/plonk/tree/kobigurk/port_to_zexe) of PLONK with arkworks backend by Kobi Gurkan


## Licensing

This code is licensed under Mozilla Public License Version 2.0 (MPL-2.0). Please see [LICENSE](https://github.com/rust-zkp/ark-plonk/blob/master/LICENSE) for further info.

## About
Initial [implementation](https://github.com/dusk-network/plonk) created by [Kev](https://github.com/kevaundray), [Carlos](https://github.com/CPerezz) and [Luke](https://github.com/LukePearson1) at Dusk Network. 
Redesigned by the [rust zkp](https://github.com/rust-zkp) team to have a backend which is compatible with the [arkworks](https://github.com/arkworks-rs) suite. This allows us to leverage the multitude of curves 
and optimised algebra present in various arkworks repositories. 

## Contributing

- If you want to contribute to this repository/project please, check [CONTRIBUTING.md](https://github.com/rust-zkp/ark-plonk/blob/master/CONTRIBUTING.md)
- If you want to report a bug or request a new feature addition, please open an issue on this repository.
