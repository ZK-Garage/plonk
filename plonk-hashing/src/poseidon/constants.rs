use crate::poseidon::{
    matrix::Matrix,
    mds::{factor_to_sparse_matrixes, MdsMatrices, SparseMatrix},
    preprocessing::compress_round_constants,
    round_constant::generate_constants,
    round_numbers::calc_round_numbers,
};
use ark_ff::PrimeField;

#[derive(Clone, Debug, PartialEq)]
pub struct PoseidonConstants<F: PrimeField> {
    pub mds_matrices: MdsMatrices<F>,
    pub round_constants: Vec<F>,
    pub compressed_round_constants: Vec<F>,
    pub pre_sparse_matrix: Matrix<F>,
    pub sparse_matrixes: Vec<SparseMatrix<F>>,
    pub domain_tag: F,
    pub full_rounds: usize,
    pub half_full_rounds: usize,
    pub partial_rounds: usize,
}

impl<F: PrimeField> PoseidonConstants<F> {
    /// Generate all constants needed for poseidon hash of specified
    /// width.  Note that WIDTH = ARITY + 1
    pub fn generate<const WIDTH: usize>() -> Self {
        let arity = WIDTH - 1;
        let mds_matrices = MdsMatrices::new(WIDTH);
        let (num_full_rounds, num_partial_rounds) =
            calc_round_numbers(WIDTH, true);

        debug_assert_eq!(num_full_rounds % 2, 0);
        let num_half_full_rounds = num_full_rounds / 2;
        let round_constants = generate_constants(
            1, // prime field
            1, // sbox
            F::size_in_bits() as u16,
            WIDTH.try_into().expect("WIDTH is too large"),
            num_full_rounds
                .try_into()
                .expect("num_full_rounds is too large"),
            num_partial_rounds
                .try_into()
                .expect("num_partial_rounds is too large"),
        );
        let domain_tag = F::from(((1 << arity) - 1) as u64);

        let compressed_round_constants = compress_round_constants(
            WIDTH,
            num_full_rounds,
            num_partial_rounds,
            &round_constants,
            &mds_matrices,
        );

        let (pre_sparse_matrix, sparse_matrixes) = factor_to_sparse_matrixes(
            mds_matrices.m.clone(),
            num_partial_rounds,
        );

        assert!(
            WIDTH * (num_full_rounds + num_partial_rounds)
                <= round_constants.len(),
            "Not enough round constants"
        );

        PoseidonConstants {
            mds_matrices,
            round_constants,
            domain_tag,
            full_rounds: num_full_rounds,
            half_full_rounds: num_half_full_rounds,
            partial_rounds: num_partial_rounds,
            compressed_round_constants,
            pre_sparse_matrix,
            sparse_matrixes,
        }
    }
}
