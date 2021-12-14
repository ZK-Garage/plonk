# ARK-PLONK library


ARK-PLONK is a generic Rust PLONK implementation using arkworks as a backend. ARK-PLONK is one of many projects implementing PLONK like: TurboPlonk, UltraPlonk. DuskPlonk, Plonky, ShPlonk, PLOOKUP, PLONKUP, halo2...etc.

ARK-PLONK is however the only generic implementation which allows any curve implementation or commitment scheme to be used and isn’t restricted to only one implementation like other existing libraries (for example duskPlonk is based on bls12-281).




## State of the art

| Library     | Description | Prover speed| Verifier speed| Proof size | Setup     |
|             |             |             |               |            |           |
| ----------- | ----------- |-------------|-------------  | -----------|-----------|
|             |             |             |               |            |           |
|             |             |             |               |            |           |
|             |             |             |               |            |           |
|             |             |             |               |            |           |
|             |             |             |               |            |           |

## Circuit implementation

ARK-PLONK's implementation is an optimization of the original PLONK protocol as it enables lookup table to the PLONK circuit. This optimization allows for precomputation of some of the operations that are not snark friendly like bit operations (see [PLOOKUP](https://eprint.iacr.org/2020/315.pdf) for further explanation on PLONK + LOOKUP tables).

Our implementation also uses custom gates similarly to [TurboPolnk](https://docs.zkproof.org/pages/standards/accepted-workshop3/proposal-turbo_plonk.pdf) which allow us to define our own custom bit arithmetic operations like efficient Poseidon or MIMC hashes which are extremely efficient to evaluate inside of a snark. 


### Proof generation
### Prover
### Verifier

### Gadgets



## Parameters

### Elliptic curve: 
Circuits in ARK-PLONK depend on two generic parameters: 


* The pairing engine which is a pairing friendly curve used for pairing operations and proof verification


* Twisted Edwards curve: derived from the first one and is used to make elliptic curve operations inside of the circuit.

  **Example:**
```rust
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
```

  Currently tests are only run with two pairing friendly curves [Bls12_381](https://lib.rs/crates/ark-bls12-381) and 
  [Bls12_377](https://docs.rs/ark-ed-on-bls12-377/0.3.0/ark_ed_on_bls12_377/#) in order to check that the library is generic and can in fact work 
  correctly with different parameters and also to measure performance when   changing the used curve.
  

         
* Commitment scheme: The first implementation of ARK-PLONK used the KZG10 commitment scheme which needs a trusted setup as explained in the KZG10 section. However, in order for other projects who don’t wish to have a trusted setup like Mithril for example there is a generic implementation using  [ark-poly-commit](https://docs.rs/ark-poly-commit/0.3.0/ark_poly_commit/) library. 



* Signature scheme: 
* Hash function:



### Simple example

In order to show how our implementation works, we will take a simple example of a gadget $a+b=c$. The gadget proves the knowledge of two private inputs $a$ and $b$ while taking $c$ as a public input.
```rust
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

```
We define three checks that our circuit relies on:
* **Range checks:** $a<26 and b<25$
* **Addition checks:** $a+b=c$
* **Elliptic curve multiplication** 
 

## Serialization


## Performance: 
To run the benchmarks and get a full report on the performance using your machine 

```RUSTFLAGS='-C target-cpu=native' cargo bench```

The benchmarks for prover and verifier are done using the [Criterion.rs](https://bheisler.github.io/criterion.rs/book/criterion_rs.html) micro-benchmarking tool. Benchmarks are repeated 10 times each and so far only been run with Bls12_381. Benchmarks taken on Intel(R) Core(TM) i9-10885H. Results are in the following table:



| Circuit size    | Prover speed | Verifier speed| 
|-----------------|--------------| --------------|            
|    $2^5$        |    9.5398ms  |    4.2881ms   |               
|    $2^6$        |    13.013ms  |    4.2781ms   |                   
|    $2^7$        |    18.137ms  |    4.2973ms   |                      
|    $2^8$        |    29.914ms  |    4.2593ms   |                       
|    $2^9$        |    50.221ms  |    4.3023ms   |                          
|    $2^{10}$     |    68.704ms  |    4.2228ms   |               
|    $2^{11}$     |    127.49ms  |    4.1379ms   |               
|    $2^{12}$     |    245.48ms  |    4.1467ms   |               
|    $2^{13}$     |    440.64ms  |    4.1770ms   |               
|    $2^{14}$     |    869.66ms  |    4.1870ms   |               
|    $2^{15}$     |    1.7712s   |    4.3390ms   |               
|    $2^{16}$     |    3.4499s   |    4.5020ms   |               
|    $2^{17}$     |    6.7577s   |    5.1572ms   |              
|    $2^{18}$     |    13.704s   |    6.8124ms   |                                     


The ark-plonk benchmarks are outperforming those of dusk-plonk library. For example, in dusk-plonk the results are taken with Intel(R) Core(TM) i9-9900X for a circuit of size $2^{16}$ are:

* Prover time: $5.46s$ which is higher than $3.44s$ in the case of ark-plonk 
* Verifier time: $9.34ms$ which is higher than $4.50ms$ in the case of Plonk.
