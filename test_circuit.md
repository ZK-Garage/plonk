```
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