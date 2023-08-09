//! acknowledgement: adapted from FileCoin Project: https://github.com/filecoin-project/neptune/blob/master/src/preprocessing.rs

use super::{matrix::vec_add, mds::MdsMatrices};
use ark_ff::vec::Vec;
use ark_ff::PrimeField;

// - Compress constants by pushing them back through linear layers and through
//   the identity components of partial layers.
// - As a result, constants need only be added after each S-box.
pub(crate) fn compress_round_constants<F: PrimeField>(
    width: usize,
    full_rounds: usize,
    partial_rounds: usize,
    round_constants: &[F],
    mds_matrices: &MdsMatrices<F>,
) -> Vec<F> {
    let inverse_matrix = &mds_matrices.m_inv;

    let mut res: Vec<F> = Vec::new();

    let round_keys = |r: usize| &round_constants[r * width..(r + 1) * width];

    // This is half full-rounds.
    let half_full_rounds = full_rounds / 2;

    // First round constants are unchanged.
    res.extend(round_keys(0));

    // Post S-box adds for the first set of full rounds should be 'inverted'
    // from next round. The final round is skipped when fully preprocessing
    // because that value must be obtained from the result of preprocessing
    // the partial rounds.
    let end = half_full_rounds - 1;
    for i in 0..end {
        let next_round = round_keys(i + 1);
        let inverted = inverse_matrix.right_apply(next_round);
        res.extend(inverted);
    }

    // The plan:
    // - Work backwards from last row in this group
    // - Invert the row.
    // - Save first constant (corresponding to the one S-box performed).
    // - Add inverted result to previous row.
    // - Repeat until all partial round key rows have been consumed.
    // - Extend the preprocessed result by the final resultant row.
    // - Move the accumulated list of single round keys to the preprocesed
    //   result.
    // - (Last produced should be first applied, so either pop until empty, or
    //   reverse and extend, etc.)

    // 'partial_keys' will accumulated the single post-S-box constant for each
    // partial-round, in reverse order.
    let mut partial_keys: Vec<F> = Vec::new();

    let final_round = half_full_rounds + partial_rounds;
    let final_round_key = round_keys(final_round).to_vec();

    // 'round_acc' holds the accumulated result of inverting and adding
    // subsequent round constants (in reverse).
    let round_acc = (0..partial_rounds)
        .map(|i| round_keys(final_round - i - 1))
        .fold(final_round_key, |acc, previous_round_keys| {
            let mut inverted = inverse_matrix.right_apply(&acc);

            partial_keys.push(inverted[0]);
            inverted[0] = F::zero();

            vec_add(previous_round_keys, &inverted)
        });

    res.extend(inverse_matrix.right_apply(&round_acc));

    while let Some(x) = partial_keys.pop() {
        res.push(x)
    }

    // Post S-box adds for the first set of full rounds should be 'inverted'
    // from next round.
    for i in 1..(half_full_rounds) {
        let start = half_full_rounds + partial_rounds;
        let next_round = round_keys(i + start);
        let inverted = inverse_matrix.right_apply(next_round);
        res.extend(inverted);
    }

    res
}
