use ark_ec::TEModelParameters;
use ark_ff::PrimeField;

use super::{StandardComposer, Variable};

pub(crate) const SBOX_ALPHA: u64 = 5;

impl<F, P> StandardComposer<F, P>
where
    F: PrimeField,
    P: TEModelParameters<BaseField = F>,
{
    /// Input 3 inputs \[w1, w2, w3\],
    /// four selectors \[q1, q2, q3, qc, qo\],
    /// constraint the following statement:
    /// ``` ignore
    ///     q1 w^5 + q2 w2^5 + q3 w3^5 + q4 + q5 w4 = 0
    /// ```
    /// Return w4
    pub fn full_affine_transform_gate(
        &mut self,
        vars: &[Variable; 3],
        selectors: &[F; 5],
    ) -> Variable {
        let w4_val = (selectors[0]
            * self.value_of_var(vars[0]).pow([SBOX_ALPHA])
            + selectors[1] * self.value_of_var(vars[1]).pow([SBOX_ALPHA])
            + selectors[2] * self.value_of_var(vars[2]).pow([SBOX_ALPHA])
            + selectors[3])
            / -selectors[4];
        let w4_var = self.add_input(w4_val);

        // add wires
        self.w_l.push(vars[0]);
        self.w_r.push(vars[1]);
        self.w_o.push(w4_var);
        self.w_4.push(vars[2]);

        // add high degree selectors
        self.q_hl.push(selectors[0]);
        self.q_hr.push(selectors[1]);
        self.q_h4.push(selectors[2]);
        self.q_c.push(selectors[3]);
        self.q_o.push(selectors[4]);

        // add dummy selector vectors
        self.q_l.push(F::zero());
        self.q_r.push(F::zero());
        self.q_m.push(F::zero());
        self.q_4.push(F::zero());
        self.q_arith.push(F::one());

        self.q_range.push(F::zero());
        self.q_logic.push(F::zero());
        self.q_fixed_group_add.push(F::zero());
        self.q_variable_group_add.push(F::zero());
        self.q_lookup.push(F::zero());

        self.perm
            .add_variables_to_map(vars[0], vars[1], w4_var, vars[2], self.n);
        self.n += 1;

        w4_var
    }

    /// Input 3 inputs \[w1, w2, w3\],
    /// four selectors \[q1, q2, q3, qc, qo\],
    /// constraint the following statement:
    /// ``` ignore
    ///     q1 w^5 + q2 w2^5 + q3 w3^5 + q4 + q5 w4 = 0
    /// ```
    /// Return w4
    pub fn partial_affine_transform_gate(
        &mut self,
        vars: &[Variable; 3],
        selectors: &[F; 5],
    ) -> Variable {
        let w4_val = (selectors[0]
            * self.value_of_var(vars[0]).pow([SBOX_ALPHA])
            + selectors[1] * self.value_of_var(vars[1])
            + selectors[2] * self.value_of_var(vars[2])
            + selectors[3])
            / -selectors[4];
        let w4_var = self.add_input(w4_val);

        // add wires
        self.w_l.push(vars[0]);
        self.w_r.push(vars[1]);
        self.w_o.push(w4_var);
        self.w_4.push(vars[2]);

        // add high degree selectors
        self.q_hl.push(selectors[0]);
        self.q_hr.push(F::zero());
        self.q_h4.push(F::zero());
        self.q_c.push(selectors[3]);
        self.q_o.push(selectors[4]);

        // add dummy selector vectors
        self.q_l.push(F::zero());
        self.q_r.push(selectors[1]);
        self.q_m.push(F::zero());
        self.q_4.push(selectors[2]);
        self.q_arith.push(F::one());

        self.q_range.push(F::zero());
        self.q_logic.push(F::zero());
        self.q_fixed_group_add.push(F::zero());
        self.q_variable_group_add.push(F::zero());
        self.q_lookup.push(F::zero());

        self.perm
            .add_variables_to_map(vars[0], vars[1], w4_var, vars[2], self.n);
        self.n += 1;

        w4_var
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        batch_test, commitment::HomomorphicCommitment,
        constraint_system::helper::gadget_tester,
    };
    use ark_bls12_377::Bls12_377;
    use ark_bls12_381::Bls12_381;

    fn test_degree_5_gates<F, P, PC>()
    where
        F: PrimeField,
        P: TEModelParameters<BaseField = F>,
        PC: HomomorphicCommitment<F>,
    {
        let res = gadget_tester::<F, P, PC>(
            |composer: &mut StandardComposer<F, P>| {
                let a = composer.add_input(F::from(2u64));
                let b = composer.add_input(F::from(3u64));
                let c = composer.add_input(F::from(5u64));

                let q1 = F::from(7u64);
                let q2 = F::from(11u64);
                let q3 = F::from(13u64);
                let q4 = F::from(17u64);
                let q5 = -F::from(1u64);

                // 7*2^5+ 11*3^5 + 13*5^5 + 17
                let d = composer.add_input(F::from(43539u64));
                // 7*2^5+ 11*3 + 13*5 + 17
                let e = composer.add_input(F::from(339u64));

                let d_rec = composer.full_affine_transform_gate(
                    &[a, b, c],
                    &[q1, q2, q3, q4, q5],
                );
                composer.assert_equal(d, d_rec);

                let e_rec = composer.partial_affine_transform_gate(
                    &[a, b, c],
                    &[q1, q2, q3, q4, q5],
                );
                composer.assert_equal(e, e_rec);

                composer.check_circuit_satisfied();
            },
            32,
        );
        assert!(res.is_ok(), "{:?}", res.err().unwrap());
    }

    // Tests for Bls12_381
    batch_test!(
        [
            test_degree_5_gates
        ],
        [] => (
            Bls12_381,
            ark_ed_on_bls12_381::EdwardsParameters
        )
    );

    // Tests for Bls12_377
    batch_test!(
        [
            test_degree_5_gates
        ],
        [] => (
            Bls12_377,
            ark_ed_on_bls12_377::EdwardsParameters
        )
    );
}
