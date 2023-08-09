//! Correct, Naive, reference implementation of Poseidon hash function.

use crate::poseidon::PoseidonError;

use crate::poseidon::constants::PoseidonConstants;
use ark_ec::TEModelParameters;
use ark_ff::PrimeField;
use core::{fmt::Debug, marker::PhantomData};
use derivative::Derivative;
use plonk_core::{constraint_system::StandardComposer, prelude as plonk};

pub trait PoseidonRefSpec<COM, const WIDTH: usize> {
    /// Field used as state
    type Field: Debug + Clone;
    /// Field used as constant paramater
    type ParameterField: PrimeField; // TODO: for now, only prime field is supported. Can be used for arkplonk
                                     // and arkworks which uses the same
                                     // PrimeField. For other field, we are not
                                     // supporting yet.

    fn full_round(
        c: &mut COM,
        constants: &PoseidonConstants<Self::ParameterField>,
        constants_offset: &mut usize,
        state: &mut [Self::Field; WIDTH],
    ) {
        let pre_round_keys = constants
            .round_constants
            .iter()
            .skip(*constants_offset)
            .map(Some);

        state.iter_mut().zip(pre_round_keys).for_each(|(l, pre)| {
            *l = Self::quintic_s_box(c, l.clone(), pre.copied(), None);
        });

        *constants_offset += WIDTH;

        Self::product_mds(c, constants, state);
    }

    fn partial_round(
        c: &mut COM,
        constants: &PoseidonConstants<Self::ParameterField>,
        constants_offset: &mut usize,
        state: &mut [Self::Field; WIDTH],
    ) {
        // TODO: we can combine add_round_constants and s_box using fewer
        // constraints
        Self::add_round_constants(c, state, constants, constants_offset);

        // apply quintic s-box to the first element
        state[0] = Self::quintic_s_box(c, state[0].clone(), None, None);

        // Multiply by MDS
        Self::product_mds(c, constants, state);
    }

    fn add_round_constants(
        c: &mut COM,
        state: &mut [Self::Field; WIDTH],
        constants: &PoseidonConstants<Self::ParameterField>,
        constants_offset: &mut usize,
    ) {
        for (element, round_constant) in state
            .iter_mut()
            .zip(constants.round_constants.iter().skip(*constants_offset))
        {
            // element.com_addi(c, round_constant);
            *element = Self::addi(c, element, round_constant)
        }

        *constants_offset += WIDTH;
    }

    fn product_mds(
        c: &mut COM,
        constants: &PoseidonConstants<Self::ParameterField>,
        state: &mut [Self::Field; WIDTH],
    ) {
        let matrix = &constants.mds_matrices.m;
        let mut result = Self::zeros::<WIDTH>(c);
        for (j, val) in result.iter_mut().enumerate() {
            for (i, row) in matrix.iter_rows().enumerate() {
                // *val += row[j] * state[i];
                let tmp = Self::muli(c, &state[i], &row[j]);
                *val = Self::add(c, val, &tmp);
            }
        }
        *state = result;
    }

    /// return (x + pre_add)^5 + post_add
    fn quintic_s_box(
        c: &mut COM,
        x: Self::Field,
        pre_add: Option<Self::ParameterField>,
        post_add: Option<Self::ParameterField>,
    ) -> Self::Field {
        let mut tmp = match pre_add {
            Some(a) => Self::addi(c, &x, &a),
            None => x.clone(),
        };
        tmp = Self::power_of_5(c, &tmp);
        match post_add {
            Some(a) => Self::addi(c, &tmp, &a),
            None => tmp,
        }
    }

    fn power_of_5(c: &mut COM, x: &Self::Field) -> Self::Field {
        let mut tmp = Self::mul(c, x, x); // x^2
        tmp = Self::mul(c, &tmp, &tmp); // x^4
        Self::mul(c, &tmp, x) // x^5
    }

    fn alloc(c: &mut COM, v: Self::ParameterField) -> Self::Field;
    fn zeros<const W: usize>(c: &mut COM) -> [Self::Field; W];
    fn zero(c: &mut COM) -> Self::Field {
        Self::zeros::<1>(c)[0].clone()
    }
    fn add(c: &mut COM, x: &Self::Field, y: &Self::Field) -> Self::Field;
    fn addi(
        c: &mut COM,
        a: &Self::Field,
        b: &Self::ParameterField,
    ) -> Self::Field;
    fn mul(c: &mut COM, x: &Self::Field, y: &Self::Field) -> Self::Field;
    fn muli(
        c: &mut COM,
        x: &Self::Field,
        y: &Self::ParameterField,
    ) -> Self::Field;
}

#[derive(Derivative)]
#[derivative(Debug(bound = ""))]
pub struct PoseidonRef<COM, S: PoseidonRefSpec<COM, WIDTH>, const WIDTH: usize>
where
    S: ?Sized,
{
    pub(crate) constants_offset: usize,
    pub(crate) current_round: usize,
    pub elements: [S::Field; WIDTH],
    pos: usize,
    pub(crate) constants: PoseidonConstants<S::ParameterField>,
}

impl<COM, S: PoseidonRefSpec<COM, WIDTH>, const WIDTH: usize>
    PoseidonRef<COM, S, WIDTH>
{
    pub fn new(
        c: &mut COM,
        constants: PoseidonConstants<S::ParameterField>,
    ) -> Self {
        let mut elements = S::zeros(c);
        elements[0] = S::alloc(c, constants.domain_tag);
        PoseidonRef {
            constants_offset: 0,
            current_round: 0,
            elements,
            pos: 1,
            constants,
        }
    }

    pub fn arity(&self) -> usize {
        WIDTH - 1
    }

    pub fn reset(&mut self, c: &mut COM) {
        self.constants_offset = 0;
        self.current_round = 0;
        self.elements[1..].iter_mut().for_each(|l| *l = S::zero(c));
        self.elements[0] = S::alloc(c, self.constants.domain_tag);
        self.pos = 1;
    }

    /// input one field element to Poseidon. Return the position of the element
    /// in state.
    pub fn input(&mut self, input: S::Field) -> Result<usize, PoseidonError> {
        // Cannot input more elements than the defined constant width
        if self.pos >= WIDTH {
            return Err(PoseidonError::FullBuffer);
        }

        // Set current element, and increase the pointer
        self.elements[self.pos] = input;
        self.pos += 1;

        Ok(self.pos - 1)
    }

    /// Output the hash
    pub fn output_hash(&mut self, c: &mut COM) -> S::Field {
        S::full_round(
            c,
            &self.constants,
            &mut self.constants_offset,
            &mut self.elements,
        );

        for _ in 1..self.constants.half_full_rounds {
            S::full_round(
                c,
                &self.constants,
                &mut self.constants_offset,
                &mut self.elements,
            );
        }

        S::partial_round(
            c,
            &self.constants,
            &mut self.constants_offset,
            &mut self.elements,
        );

        for _ in 1..self.constants.partial_rounds {
            S::partial_round(
                c,
                &self.constants,
                &mut self.constants_offset,
                &mut self.elements,
            );
        }

        for _ in 0..self.constants.half_full_rounds {
            S::full_round(
                c,
                &self.constants,
                &mut self.constants_offset,
                &mut self.elements,
            )
        }

        self.elements[1].clone()
    }
}

pub struct NativeSpecRef<F: PrimeField> {
    _field: PhantomData<F>,
}

impl<F: PrimeField, const WIDTH: usize> PoseidonRefSpec<(), WIDTH>
    for NativeSpecRef<F>
{
    type Field = F;
    type ParameterField = F;

    fn alloc(_c: &mut (), v: Self::ParameterField) -> Self::Field {
        v
    }

    fn zeros<const W: usize>(_c: &mut ()) -> [Self::Field; W] {
        [F::zero(); W]
    }

    fn add(_c: &mut (), x: &Self::Field, y: &Self::Field) -> Self::Field {
        *x + *y
    }

    fn addi(
        _c: &mut (),
        a: &Self::Field,
        b: &Self::ParameterField,
    ) -> Self::Field {
        *a + *b
    }

    fn mul(_c: &mut (), x: &Self::Field, y: &Self::Field) -> Self::Field {
        *x * *y
    }

    fn muli(
        _c: &mut (),
        x: &Self::Field,
        y: &Self::ParameterField,
    ) -> Self::Field {
        *x * *y
    }
}

pub struct PlonkSpecRef;

impl<F, P, const WIDTH: usize>
    PoseidonRefSpec<plonk::StandardComposer<F, P>, WIDTH> for PlonkSpecRef
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    type Field = plonk::Variable;
    type ParameterField = F;

    fn alloc(
        c: &mut StandardComposer<F, P>,
        v: Self::ParameterField,
    ) -> Self::Field {
        c.add_input(v)
    }

    fn zeros<const W: usize>(
        c: &mut StandardComposer<F, P>,
    ) -> [Self::Field; W] {
        [c.zero_var(); W]
    }

    fn add(
        c: &mut StandardComposer<F, P>,
        x: &Self::Field,
        y: &Self::Field,
    ) -> Self::Field {
        c.arithmetic_gate(|g| g.witness(*x, *y, None).add(F::one(), F::one()))
    }

    fn addi(
        c: &mut StandardComposer<F, P>,
        a: &Self::Field,
        b: &Self::ParameterField,
    ) -> Self::Field {
        let zero = c.zero_var();
        c.arithmetic_gate(|g| {
            g.witness(*a, zero, None)
                .add(F::one(), F::zero())
                .constant(*b)
        })
    }

    fn mul(
        c: &mut StandardComposer<F, P>,
        x: &Self::Field,
        y: &Self::Field,
    ) -> Self::Field {
        c.arithmetic_gate(|q| q.witness(*x, *y, None).mul(F::one()))
    }

    fn muli(
        c: &mut StandardComposer<F, P>,
        x: &Self::Field,
        y: &Self::ParameterField,
    ) -> Self::Field {
        let zero = c.zero_var();
        c.arithmetic_gate(|g| g.witness(*x, zero, None).add(*y, F::zero()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::PairingEngine;

    type E = ark_bls12_381::Bls12_381;
    type P = ark_ed_on_bls12_381::EdwardsParameters;
    type Fr = <E as PairingEngine>::Fr;
    use ark_std::{test_rng, UniformRand};

    #[test]
    // poseidon should output something if num_inputs = arity
    fn sanity_test() {
        const ARITY: usize = 4;
        const WIDTH: usize = ARITY + 1;
        let mut rng = test_rng();

        let param = PoseidonConstants::generate::<WIDTH>();
        let mut poseidon = PoseidonRef::<(), NativeSpecRef<Fr>, WIDTH>::new(
            &mut (),
            param.clone(),
        );
        let inputs = (0..ARITY).map(|_| Fr::rand(&mut rng)).collect::<Vec<_>>();

        inputs.iter().for_each(|x| {
            let _ = poseidon.input(*x).unwrap();
        });
        let native_hash: Fr = poseidon.output_hash(&mut ());

        let mut c = StandardComposer::<Fr, P>::new();
        let inputs_var =
            inputs.iter().map(|x| c.add_input(*x)).collect::<Vec<_>>();
        let mut poseidon_circuit =
            PoseidonRef::<_, PlonkSpecRef, WIDTH>::new(&mut c, param);
        inputs_var.iter().for_each(|x| {
            let _ = poseidon_circuit.input(*x).unwrap();
        });
        let plonk_hash = poseidon_circuit.output_hash(&mut c);

        c.check_circuit_satisfied();

        let expected = c.add_input(native_hash);
        c.assert_equal(expected, plonk_hash);

        c.check_circuit_satisfied();
        println!(
            "circuit size for WIDTH {} poseidon: {}",
            WIDTH,
            c.circuit_bound()
        )
    }

    #[test]
    #[should_panic]
    // poseidon should output something if num_inputs > arity
    fn sanity_test_failure() {
        const ARITY: usize = 4;
        const WIDTH: usize = ARITY + 1;
        let mut rng = test_rng();

        let param = PoseidonConstants::generate::<WIDTH>();
        let mut poseidon =
            PoseidonRef::<(), NativeSpecRef<Fr>, WIDTH>::new(&mut (), param);
        (0..(ARITY + 1)).for_each(|_| {
            let _ = poseidon.input(Fr::rand(&mut rng)).unwrap();
        });
        let _ = poseidon.output_hash(&mut ());
    }

    // use crate::tests::{
    //     conversion::cast_field,
    //     neptune_hyper_parameter::collect_neptune_constants,
    // };
    // use neptune::{
    //     poseidon::{HashMode, PoseidonConstants as NeptunePoseidonConstants},
    //     Strength,
    // };
    // let constants = NeptunePoseidonConstants::<Fr,
    // A>::new_with_strength(strength); let mut p = NeptunePoseidon::<Fr,
    // A>::new(&constants); let mut p2 = NeptunePoseidon::<Fr,
    // A>::new(&constants); let mut p3 = NeptunePoseidon::<Fr,
    // A>::new(&constants); let mut p4 = NeptunePoseidon::<Fr,
    // A>::new(&constants);

    //     let test_arity = constants.arity();
    //     for n in 0..test_arity {
    //         let scalar = Fr::from(n as u64);
    //         p.input(scalar).unwrap();
    //         p2.input(scalar).unwrap();
    //         p3.input(scalar).unwrap();
    //         p4.input(scalar).unwrap();
    //     }

    //     let digest = p.hash();
    //     let digest2 = p2.hash_in_mode(Correct);
    //     let digest3 = p3.hash_in_mode(OptimizedStatic);
    //     let digest4 = p4.hash_in_mode(OptimizedDynamic);

    // #[test]
    // fn compare_with_neptune() {
    //     const ARITY: usize = 2;
    //     const WIDTH: usize = ARITY + 1;
    //     type NepArity = generic_array::typenum::U2;

    //     let (nep_consts, ark_consts) =
    //         collect_neptune_constants::<NepArity>(Strength::Standard);

    //     let mut rng = test_rng();
    //     let inputs_ff = (0..ARITY)
    //         .map(|_| blstrs::Scalar::random(&mut rng))
    //         .collect::<Vec<_>>();
    //     let inputs =
    //         inputs_ff.iter().map(|&x| cast_field(x)).collect::<Vec<_>>();

    //     let mut neptune_poseidon =
    //         neptune::Poseidon::<blstrs::Scalar, NepArity>::new(&nep_consts);
    //     let mut ark_poseidon = PoseidonRef::<(), NativeSpecRef<Fr>,
    // WIDTH>::new(         &mut (),
    //         ark_consts,
    //     );

    //     inputs_ff.iter().for_each(|x| {
    //         neptune_poseidon.input(*x).unwrap();
    //     });
    //     inputs.iter().for_each(|x| {
    //         ark_poseidon.input(*x).unwrap();
    //     });

    //     let digest_expected =
    //         cast_field(neptune_poseidon.hash_in_mode(HashMode::Correct));
    //     let digest_actual = ark_poseidon.output_hash(&mut ());

    //     assert_eq!(digest_expected, digest_actual);
    // }
}
