// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! Constants are used in the permutation argument to generate H cosets.
#![allow(non_snake_case)]

use ark_ff::FftField;

pub(crate) fn K1<F: FftField>() -> F {
    F::from(7_u64)
}

pub(crate) fn K2<F: FftField>() -> F {
    F::from(13_u64)
}

pub(crate) fn K3<F: FftField>() -> F {
    F::from(17_u64)
}

#[cfg(test)]
mod test {
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;
    use ark_poly::{EvaluationDomain, GeneralEvaluationDomain};

    use crate::batch_test_field;

    use super::*;

    /// Check if `cts` generate valid cosets of the roots of
    /// unity subgroup (of `domain_size`) of the field F.
    /// https://hackmd.io/CfFCbA0TTJ6X08vHg0-9_g
    fn check_constants<F>(cts: &[F], domain_size: u64) -> bool
    where
        F: FftField,
    {
        if cts.is_empty() {
            return true;
        };

        let domain: GeneralEvaluationDomain<F> =
            EvaluationDomain::new(domain_size as usize).unwrap();

        // First check
        if domain.evaluate_vanishing_polynomial(cts[0]).is_zero() {
            return false;
        }

        let mut prev_cts = Vec::with_capacity(cts.len());
        prev_cts.push(cts[0]);

        // Rest of the constants
        for k_last in cts.iter().skip(1) {
            let k_last_inv = k_last.inverse().unwrap();

            // Check that the constant k_last is not in each of the previously
            // generated cosets k_prev * H.
            // <=> (k_prev / k_last )^domain_size  -1 != 0.
            if prev_cts.iter().any(|&k_prev| {
                domain
                    .evaluate_vanishing_polynomial(k_prev * k_last_inv)
                    .is_zero()
            }) {
                return false;
            }
            prev_cts.push(*k_last);
        }

        true
    }

    fn test_constants<F>()
    where
        F: FftField,
    {
        // The constants are checked for subgropus H up to 2^MAX_DEGREE size.
        let MAX_DEGREE = 32;
        let constants: [F; 3] = [super::K1(), super::K2(), super::K3()];
        assert!(check_constants(&constants, 2u64.pow(MAX_DEGREE)));
    }

    batch_test_field!(
        [test_constants],
        [] => (
            Bls12_381
        )
    );
    batch_test_field!(
        [test_constants],
        [] => (
            Bls12_377
        )
    );
}
