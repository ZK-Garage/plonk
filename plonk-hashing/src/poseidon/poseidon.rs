//! optimized poseidon

use crate::poseidon::constants::PoseidonConstants;
use crate::poseidon::matrix::Matrix;
use crate::poseidon::mds::SparseMatrix;
use crate::poseidon::PoseidonError;
use ark_ec::TEModelParameters;
use ark_ff::PrimeField;
use core::{fmt::Debug, marker::PhantomData};
use derivative::Derivative;
use plonk_core::constraint_system::{StandardComposer, Variable};
use plonk_core::prelude as plonk;

// TODO: reduce duplicate code with `poseidon_ref`
pub trait PoseidonSpec<COM, const WIDTH: usize> {
    type Field: Debug + Clone;
    type ParameterField: PrimeField;

    fn output_hash(
        c: &mut COM,
        elements: &[Self::Field; WIDTH],
        constants: &PoseidonConstants<Self::ParameterField>,
    ) -> Self::Field {
        // Counters
        let mut constants_offset = 0usize;
        let mut current_round = 0usize;
        // State vector to modify
        let mut state = elements.clone();
        Self::add_round_constants(
            c,
            &mut state,
            constants,
            &mut constants_offset,
        );

        for _ in 0..constants.half_full_rounds {
            Self::full_round(
                c,
                constants,
                &mut current_round,
                &mut constants_offset,
                false,
                &mut state,
            )
        }

        for _ in 0..constants.partial_rounds {
            Self::partial_round(
                c,
                constants,
                &mut current_round,
                &mut constants_offset,
                &mut state,
            );
        }

        // All but last full round
        for _ in 1..constants.half_full_rounds {
            Self::full_round(
                c,
                constants,
                &mut current_round,
                &mut constants_offset,
                false,
                &mut state,
            );
        }
        Self::full_round(
            c,
            constants,
            &mut current_round,
            &mut constants_offset,
            true,
            &mut state,
        );

        assert_eq!(
            constants_offset,
            constants.compressed_round_constants.len(),
            "Constants consumed ({}) must equal preprocessed constants provided ({}).",
            constants_offset,
            constants.compressed_round_constants.len()
        );

        state[1].clone()
    }

    fn full_round(
        c: &mut COM,
        constants: &PoseidonConstants<Self::ParameterField>,
        current_round: &mut usize,
        const_offset: &mut usize,
        last_round: bool,
        state: &mut [Self::Field; WIDTH],
    ) {
        let to_take = WIDTH;
        let post_round_keys = constants
            .compressed_round_constants
            .iter()
            .skip(*const_offset)
            .take(to_take);

        if !last_round {
            let needed = *const_offset + to_take;
            assert!(
                needed <= constants.compressed_round_constants.len(),
                "Not enough preprocessed round constants ({}), need {}.",
                constants.compressed_round_constants.len(),
                needed
            );
        }

        state.iter_mut().zip(post_round_keys).for_each(|(l, post)| {
            // Be explicit that no round key is added after last round of S-boxes.
            let post_key = if last_round {
                panic!(
                    "Trying to skip last full round, but there is a key here! ({:?})",
                    post
                );
            } else {
                Some(post.clone())
            };
            *l = Self::quintic_s_box(c, l.clone(), None, post_key);
        });

        if last_round {
            state.iter_mut().for_each(|l| {
                *l = Self::quintic_s_box(c, l.clone(), None, None)
            })
        } else {
            *const_offset += to_take;
        }
        Self::round_product_mds(c, constants, current_round, state);
    }

    fn partial_round(
        c: &mut COM,
        constants: &PoseidonConstants<Self::ParameterField>,
        current_round: &mut usize,
        const_offset: &mut usize,
        state: &mut [Self::Field; WIDTH],
    ) {
        let post_round_key =
            constants.compressed_round_constants[*const_offset];

        state[0] = Self::quintic_s_box(
            c,
            state[0].clone(),
            None,
            Some(post_round_key),
        );
        *const_offset += 1;

        Self::round_product_mds(c, constants, current_round, state);
    }

    fn add_round_constants(
        c: &mut COM,
        state: &mut [Self::Field; WIDTH],
        constants: &PoseidonConstants<Self::ParameterField>,
        const_offset: &mut usize,
    ) {
        for (element, round_constant) in state.iter_mut().zip(
            constants
                .compressed_round_constants
                .iter()
                .skip(*const_offset),
        ) {
            *element = Self::addi(c, element, round_constant);
        }
        *const_offset += WIDTH;
    }

    fn round_product_mds(
        c: &mut COM,
        constants: &PoseidonConstants<Self::ParameterField>,
        current_round: &mut usize,
        state: &mut [Self::Field; WIDTH],
    ) {
        let full_half = constants.half_full_rounds;
        let sparse_offset = full_half - 1;
        if *current_round == sparse_offset {
            Self::product_mds_with_matrix(
                c,
                state,
                &constants.pre_sparse_matrix,
            )
        } else {
            if (*current_round > sparse_offset)
                && (*current_round < full_half + constants.partial_rounds)
            {
                let index = *current_round - sparse_offset - 1;
                let sparse_matrix = &constants.sparse_matrixes[index];

                Self::product_mds_with_sparse_matrix(c, state, sparse_matrix)
            } else {
                Self::product_mds(c, constants, state)
            }
        };

        *current_round += 1;
    }

    fn product_mds(
        c: &mut COM,
        constants: &PoseidonConstants<Self::ParameterField>,
        state: &mut [Self::Field; WIDTH],
    ) {
        Self::product_mds_with_matrix(c, state, &constants.mds_matrices.m)
    }

    fn linear_combination(
        c: &mut COM,
        state: &[Self::Field; WIDTH],
        coeff: impl IntoIterator<Item = Self::ParameterField>,
    ) -> Self::Field {
        state.iter().zip(coeff).fold(Self::zero(c), |acc, (x, y)| {
            let tmp = Self::muli(c, x, &y);
            Self::add(c, &tmp, &acc)
        })
    }

    /// compute state @ Mat where `state` is a row vector
    fn product_mds_with_matrix(
        c: &mut COM,
        state: &mut [Self::Field; WIDTH],
        matrix: &Matrix<Self::ParameterField>,
    ) {
        let mut result = Self::zeros::<WIDTH>(c);
        for (col_index, val) in result.iter_mut().enumerate() {
            // for (i, row) in matrix.iter_rows().enumerate() {
            //     // *val += row[j] * state[i];
            //     let tmp = Self::muli(c, &state[i], &row[j]);
            //     *val = Self::add(c, val, &tmp);
            // }
            *val = Self::linear_combination(
                c,
                state,
                matrix.column(col_index).cloned(),
            );
        }

        *state = result;
    }

    fn product_mds_with_sparse_matrix(
        c: &mut COM,
        state: &mut [Self::Field; WIDTH],
        matrix: &SparseMatrix<Self::ParameterField>,
    ) {
        let mut result = Self::zeros::<WIDTH>(c);

        // First column is dense.
        // for (i, val) in matrix.w_hat.iter().enumerate() {
        //     // result[0] += w_hat[i] * state[i];
        //     let tmp = Self::muli(c, &state[i], &val);
        //     result[0] = Self::add(c, &result[0], &tmp);
        // }
        result[0] =
            Self::linear_combination(c, state, matrix.w_hat.iter().cloned());

        for (j, val) in result.iter_mut().enumerate().skip(1) {
            // for each j, result[j] = state[j] + state[0] * v_rest[j-1]

            // Except for first row/column, diagonals are one.
            *val = Self::add(c, val, &state[j]);
            //
            // // First row is dense.
            let tmp = Self::muli(c, &state[0], &matrix.v_rest[j - 1]);
            *val = Self::add(c, val, &tmp);
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
pub struct Poseidon<COM, S: PoseidonSpec<COM, WIDTH>, const WIDTH: usize>
where
    S: ?Sized,
{
    pub(crate) constants: PoseidonConstants<S::ParameterField>,
}

impl<COM, S: PoseidonSpec<COM, WIDTH>, const WIDTH: usize>
    Poseidon<COM, S, WIDTH>
where
    S: ?Sized,
{
    pub fn new(constants: PoseidonConstants<S::ParameterField>) -> Self {
        Poseidon { constants }
    }

    pub fn arity(&self) -> usize {
        WIDTH - 1
    }

    /// Hash an array of ARITY-many elements.  The size of elements could be
    /// specified as WIDTH - 1 when const generic expressions are allowed.
    /// Function will panic if elements does not have length ARITY.
    pub fn output_hash(&self, elements: &[S::Field], c: &mut COM) -> S::Field {
        // Inputs should have domain_tag as its leading entry
        let mut inputs = S::zeros(c);
        // clone_from_slice will panic unless we provided ARITY-many elements
        inputs[1..WIDTH].clone_from_slice(&elements[..(WIDTH - 1)]);

        S::output_hash(c, &inputs, &self.constants)
    }
}

pub struct NativeSpec<F: PrimeField, const WIDTH: usize> {
    _field: PhantomData<F>,
}

impl<F: PrimeField, const WIDTH: usize> PoseidonSpec<(), WIDTH>
    for NativeSpec<F, WIDTH>
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

pub struct PlonkSpec<const WIDTH: usize>;

impl<F, P, const WIDTH: usize>
    PoseidonSpec<plonk::StandardComposer<F, P>, WIDTH> for PlonkSpec<WIDTH>
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

    #[cfg(not(feature = "no-optimize"))]
    fn quintic_s_box(
        c: &mut StandardComposer<F, P>,
        x: Self::Field,
        pre_add: Option<Self::ParameterField>,
        post_add: Option<Self::ParameterField>,
    ) -> Self::Field {
        match (pre_add, post_add) {
            (None, None) => Self::power_of_5(c, &x),
            (Some(_), None) => {
                unreachable!("currently no one is using this")
            }
            (None, Some(post_add)) => {
                let x_2 = Self::mul(c, &x, &x);
                let x_4 = Self::mul(c, &x_2, &x_2);
                c.arithmetic_gate(|g| {
                    g.witness(x_4, x, None).mul(F::one()).constant(post_add)
                })
            }
            (Some(_), Some(_)) => {
                unreachable!("currently no one is using this")
            }
        }
    }

    #[cfg(not(feature = "no-optimize"))]
    fn linear_combination(
        c: &mut StandardComposer<F, P>,
        state: &[Self::Field; WIDTH],
        coeff: impl IntoIterator<Item = Self::ParameterField>,
    ) -> Self::Field {
        let coeffs = coeff.into_iter().collect::<Vec<_>>();
        let mut remaining = WIDTH;
        let mut index = 0;
        let mut result: Self::Field;
        // the first time you have no accumulated result yet, so you can take 3
        // inputs
        if remaining < 3 {
            // this is unlikely, WIDTH is usually at least 3
            result = c.arithmetic_gate(|g| {
                g.witness(state[0], state[1], None)
                    .add(coeffs[0], coeffs[1])
            });
            remaining -= 2;
        } else {
            result = c.arithmetic_gate(|g| {
                g.witness(state[index], state[index + 1], None)
                    .add(coeffs[index], coeffs[index + 1])
                    .fan_in_3(coeffs[index + 2], state[index + 2])
            });
            index += 3;
            remaining -= 3;
        }

        // Now you have an accumulated result to carry, so can only take 2
        // inputs at a time
        while remaining > 0 {
            if remaining < 2 {
                // Accumulate remaining one
                result = c.arithmetic_gate(|g| {
                    g.witness(state[index], result, None)
                        .add(coeffs[index], Self::ParameterField::one())
                });
                remaining -= 1;
            } else {
                // Accumulate next two
                result = c.arithmetic_gate(|g| {
                    g.witness(state[index], state[index + 1], None)
                        .add(coeffs[index], coeffs[index + 1])
                        .fan_in_3(Self::ParameterField::one(), result)
                });
                index += 2;
                remaining -= 2;
            }
        }
        result
    }

    #[cfg(not(feature = "no-optimize"))]
    fn product_mds_with_sparse_matrix(
        c: &mut StandardComposer<F, P>,
        state: &mut [Self::Field; WIDTH],
        matrix: &SparseMatrix<Self::ParameterField>,
    ) {
        let mut result = Self::zeros::<WIDTH>(c);

        result[0] =
            Self::linear_combination(c, state, matrix.w_hat.iter().cloned());
        for (j, val) in result.iter_mut().enumerate().skip(1) {
            // for each j, result[j] = state[j] + state[0] * v_rest[j-1]
            *val = c.arithmetic_gate(|g| {
                g.witness(state[0], state[j], None)
                    .add(matrix.v_rest[j - 1], F::one())
            });
        }
        *state = result;
    }
}

trait HashFunction<const WIDTH: usize, COM = ()> {
    type Input;
    type Output;

    // The input ought to have size ARITY, but const generic expressions aren't
    // allowed yet
    fn hash(&self, input: &[Self::Input], compiler: &mut COM) -> Self::Output;
}

impl<const WIDTH: usize, F, P> HashFunction<WIDTH, StandardComposer<F, P>>
    for Poseidon<StandardComposer<F, P>, PlonkSpec<WIDTH>, WIDTH>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    type Input = Variable;
    type Output = Variable;

    fn hash(
        &self,
        input: &[Self::Input],
        compiler: &mut StandardComposer<F, P>,
    ) -> Self::Output {
        self.output_hash(input, compiler)
    }
}
