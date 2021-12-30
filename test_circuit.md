Trait that should be implemented for any circuit function to provide to it
the capabilities of automatically being able to generate, and verify proofs
as well as compile the circuit.
## Example
```rust
use ark_bls12_381::{Bls12_381, Fr as BlsScalar};
use ark_ec::PairingEngine;
use ark_ec::models::twisted_edwards_extended::GroupAffine;
use ark_ec::{TEModelParameters, AffineCurve, ProjectiveCurve};
use ark_ed_on_bls12_381::{
    EdwardsAffine as JubJubAffine, EdwardsParameters as JubJubParameters,
    EdwardsProjective as JubJubProjective, Fr as JubJubScalar,
};
use ark_ff::{PrimeField, FftField, BigInteger};
use ark_plonk::prelude::*;
use ark_poly::polynomial::univariate::DensePolynomial;
use ark_poly_commit::PolynomialCommitment;
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize};
use num_traits::{Zero, One};
use rand::rngs::OsRng;
fn main() -> Result<(), Error> {
// Implements a circuit that checks:
// 1) a + b = c where C is a PI
// 2) a <= 2^6
// 3) b <= 2^5
// 4) a * b = d where D is a PI
// 5) JubJub::GENERATOR * e(JubJubScalar) = f where F is a PI
#[derive(derivative::Derivative)]
#[derivative(Debug(bound = ""), Default(bound = ""))]
pub struct TestCircuit<
    F: FftField + PrimeField,
    P: TEModelParameters<BaseField = F>,
> {
    a: F,
    b: F,
    c: F,
    d: F,
    e: P::ScalarField,
    f: GroupAffine<P>,
}
impl<F, P> Circuit<F, P> for TestCircuit<F, P>
where
    F: FftField + PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    const CIRCUIT_ID: [u8; 32] = [0xff; 32];
    fn gadget(
        &mut self,
        composer: &mut StandardComposer<F, P>,
    ) -> Result<(), Error> {
        let a = composer.add_input(self.a);
        let b = composer.add_input(self.b);
        let zero = composer.zero_var();
        // Make first constraint a + b = c (as public input)
        composer.arithmetic_gate(|gate| {
            gate.witness(a, b, Some(zero))
                .add(F::one(), F::one())
                .pi(-self.c)
        });
        // Check that a and b are in range
        composer.range_gate(a, 1 << 6);
        composer.range_gate(b, 1 << 5);
        // Make second constraint a * b = d
        composer.arithmetic_gate(|gate| {
            gate.witness(a, b, Some(zero)).mul(F::one()).pi(-self.d)
        });
        let e = composer
            .add_input(from_embedded_curve_scalar::<F, P>(self.e));
        let (x, y) = P::AFFINE_GENERATOR_COEFFS;
        let generator = GroupAffine::new(x, y);
        let scalar_mul_result =
            composer.fixed_base_scalar_mul(e, generator);
        // Apply the constrain
        composer.assert_equal_public_point(scalar_mul_result, self.f);
        Ok(())
    }
    fn padded_circuit_size(&self) -> usize {
        1 << 11
    }
}
// Generate CRS
type PC = ark_plonk::commitment::KZG10::<Bls12_381>;
let pp = PC::setup(
    1 << 12,
    None,
    &mut OsRng,
)?;
let mut circuit = TestCircuit::<BlsScalar, JubJubParameters>::default();
// Compile the circuit
let (pk_p, verifier_data) = circuit.compile(&pp)?;
let (x, y) = JubJubParameters::AFFINE_GENERATOR_COEFFS;
let generator: GroupAffine<JubJubParameters> = GroupAffine::new(x, y);
let point_f_pi: GroupAffine<JubJubParameters> = AffineCurve::mul(
    &generator,
    JubJubScalar::from(2u64).into_repr(),
)
.into_affine();
// Prover POV
let proof = {
    let mut circuit: TestCircuit<BlsScalar, JubJubParameters> = TestCircuit {
        a: BlsScalar::from(20u64),
        b: BlsScalar::from(5u64),
        c: BlsScalar::from(25u64),
        d: BlsScalar::from(100u64),
        e: JubJubScalar::from(2u64),
        f: point_f_pi,
    };
    circuit.gen_proof(&pp, pk_p, b"Test")?
};

// Verifier POV
let public_inputs: Vec<PublicInputValue<BlsScalar>> = vec![
    BlsScalar::from(25u64).into(),
    BlsScalar::from(100u64).into(),
    GeIntoPubInput::into_pi(point_f_pi),
];
let VerifierData { key, pi_pos } = verifier_data;

verify_proof::<BlsScalar, JubJubParameters, PC>(
    &pp,
    key,
    &proof,
    &public_inputs,
    &pi_pos,
    b"Test",
)
}
```
