//! acknowledgement: adapted from FileCoin Project: https://github.com/filecoin-project/neptune/blob/master/src/matrix.rs

use ark_ff::PrimeField;
use core::ops::{Index, IndexMut};

#[derive(Clone, Eq, PartialEq, Debug)]
pub struct Matrix<T: Clone>(pub Vec<Vec<T>>);

impl<T: Clone> From<Vec<Vec<T>>> for Matrix<T> {
    fn from(v: Vec<Vec<T>>) -> Self {
        Matrix(v)
    }
}

impl<T: Clone> Matrix<T> {
    pub fn num_rows(&self) -> usize {
        self.0.len()
    }

    pub fn num_columns(&self) -> usize {
        if self.0.is_empty() {
            0
        } else {
            let column_length = self.0[0].len();
            for row in &self.0 {
                if row.len() != column_length {
                    panic!("not a matrix");
                }
            }
            column_length
        }
    }

    pub fn iter_rows(&self) -> impl Iterator<Item = &'_ Vec<T>> {
        self.0.iter()
    }

    pub fn column(&self, column: usize) -> impl Iterator<Item = &'_ T> {
        self.0.iter().map(move |row| &row[column])
    }

    pub fn is_square(&self) -> bool {
        self.num_rows() == self.num_columns()
    }

    pub fn transpose(&self) -> Matrix<T> {
        let size = self.num_rows();
        let mut new = Vec::with_capacity(size);
        for j in 0..size {
            let mut row = Vec::with_capacity(size);
            for i in 0..size {
                row.push(self.0[i][j].clone())
            }
            new.push(row);
        }
        Matrix(new)
    }
}

impl<T: Clone> Index<usize> for Matrix<T> {
    type Output = Vec<T>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T: Clone> IndexMut<usize> for Matrix<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

// from iterator rows
impl<F: Clone> FromIterator<Vec<F>> for Matrix<F> {
    fn from_iter<T: IntoIterator<Item = Vec<F>>>(iter: T) -> Self {
        let rows = iter.into_iter().collect::<Vec<_>>();
        Self(rows)
    }
}

impl<F: PrimeField> Matrix<F> {
    /// return an identity matrix of size `n*n`
    pub fn identity(n: usize) -> Matrix<F> {
        let mut m = Matrix(vec![vec![F::zero(); n]; n]);
        for i in 0..n {
            m.0[i][i] = F::one();
        }
        m
    }

    pub fn is_identity(&self) -> bool {
        if !self.is_square() {
            return false;
        }
        for i in 0..self.num_rows() {
            for j in 0..self.num_columns() {
                if self.0[i][j] != kronecker_delta(i, j) {
                    return false;
                }
            }
        }
        true
    }

    /// check if `self` is square and `self[1..][1..]` is identity
    pub fn is_sparse(&self) -> bool {
        self.is_square() && self.minor(0, 0).is_identity()
    }

    pub fn mul_by_scalar(&self, scalar: F) -> Self {
        let res = self
            .0
            .iter()
            .map(|row| {
                row.iter()
                    .map(|val| {
                        let mut prod = scalar;
                        prod.mul_assign(val);
                        prod
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        Matrix(res)
    }

    /// return `self @ vec`, treating `vec` as a column vector.
    pub fn mul_col_vec(&self, v: &[F]) -> Vec<F> {
        assert!(
            self.is_square(),
            "Only square matrix can be applied to vector."
        );
        assert_eq!(
            self.num_rows(),
            v.len(),
            "Matrix can only be applied to vector of same size."
        );

        let mut result = vec![F::zero(); v.len()];

        for (result, row) in result.iter_mut().zip(self.0.iter()) {
            for (mat_val, vec_val) in row.iter().zip(v) {
                let mut tmp = *mat_val;
                tmp.mul_assign(vec_val);
                result.add_assign(&tmp);
            }
        }
        result
    }

    /// return `vec @ self`, treat `vec` as a row vector.
    pub fn right_apply(&self, v: &[F]) -> Vec<F> {
        self.mul_row_vec_at_left(v)
    }

    /// return `self @ vec`, treating `vec` as a column vector.
    pub fn left_apply(&self, v: &[F]) -> Vec<F> {
        self.mul_col_vec(v)
    }

    /// return `vec @ self`, treating `vec` as a row vector.
    pub fn mul_row_vec_at_left(&self, v: &[F]) -> Vec<F> {
        assert!(
            self.is_square(),
            "Only square matrix can be applied to vector."
        );
        assert_eq!(
            self.num_rows(),
            v.len(),
            "Matrix can only be applied to vector of same size."
        );

        let mut result = vec![F::zero(); v.len()];
        for (j, val) in result.iter_mut().enumerate() {
            for (i, row) in self.0.iter().enumerate() {
                let mut tmp = row[j];
                tmp.mul_assign(&v[i]);
                val.add_assign(&tmp);
            }
        }
        result
    }

    /// return `self @ other`
    pub fn matmul(&self, other: &Self) -> Option<Self> {
        if self.num_rows() != other.num_columns() {
            return None;
        };

        let other_t = other.transpose();

        let res = self
            .0
            .iter()
            .map(|input_row| {
                other_t
                    .iter_rows()
                    .map(|transposed_column| {
                        inner_product(input_row, transposed_column)
                    })
                    .collect()
            })
            .collect();
        Some(Matrix(res))
    }

    pub fn invert(&self) -> Option<Self> {
        let mut shadow = Self::identity(self.num_columns());
        let ut = self.upper_triangular(&mut shadow);

        ut.and_then(|x| x.reduce_to_identity(&mut shadow))
            .and(Some(shadow))
    }

    pub fn is_invertible(&self) -> bool {
        self.is_square() && self.invert().is_some()
    }

    pub fn minor(&self, i: usize, j: usize) -> Self {
        assert!(self.is_square());
        let size = self.num_rows();
        assert!(size > 0);
        let new: Vec<Vec<F>> = self
            .0
            .iter()
            .enumerate()
            .filter_map(|(ii, row)| {
                if ii == i {
                    None
                } else {
                    let mut new_row = row.clone();
                    new_row.remove(j);
                    Some(new_row)
                }
            })
            .collect();
        let res = Matrix(new);
        assert!(res.is_square());
        res
    }

    /// Assumes matrix is partially reduced to upper triangular. `column` is the
    /// column to eliminate from all rows. Returns `None` if either:
    ///   - no non-zero pivot can be found for `column`
    ///   - `column` is not the first
    pub fn eliminate(&self, column: usize, shadow: &mut Self) -> Option<Self> {
        let zero = F::zero();
        let pivot_index = (0..self.num_rows()).find(|&i| {
            self[i][column] != zero && (0..column).all(|j| self[i][j] == zero)
        })?;

        let pivot = &self[pivot_index];
        let pivot_val = pivot[column];

        // This should never fail since we have a non-zero `pivot_val` if we got
        // here.
        let inv_pivot = pivot_val.inverse()?;
        let mut result = Vec::with_capacity(self.num_rows());
        result.push(pivot.clone());

        for (i, row) in self.iter_rows().enumerate() {
            if i == pivot_index {
                continue;
            };

            let val = row[column];
            if val == zero {
                result.push(row.to_vec());
            } else {
                let factor = val * inv_pivot;
                let scaled_pivot = scalar_vec_mul(factor, pivot);
                let eliminated = vec_sub(row, &scaled_pivot);
                result.push(eliminated);

                let shadow_pivot = &shadow[pivot_index];
                let scaled_shadow_pivot = scalar_vec_mul(factor, shadow_pivot);
                let shadow_row = &shadow[i];
                shadow[i] = vec_sub(shadow_row, &scaled_shadow_pivot);
            }
        }

        let pivot_row = shadow.0.remove(pivot_index);
        shadow.0.insert(0, pivot_row);

        Some(result.into())
    }

    /// Performs row operations to put a matrix in upper triangular form.
    /// Each row operation is performed on `shadow` as well to keep track
    /// of their cumulative effect.  In other words, row operations are
    /// performed on the augmented matrix [self | shadow ].
    pub fn upper_triangular(&self, shadow: &mut Self) -> Option<Self> {
        assert!(self.is_square());
        let mut result = Vec::with_capacity(self.num_rows());
        let mut shadow_result = Vec::with_capacity(self.num_rows());

        let mut curr = self.clone();
        let mut column = 0;
        while curr.num_rows() > 1 {
            let initial_rows = curr.num_rows();

            curr = curr.eliminate(column, shadow)?;
            result.push(curr[0].clone());
            shadow_result.push(shadow[0].clone());
            column += 1;

            curr = Matrix(curr.0[1..].to_vec());
            *shadow = Matrix(shadow.0[1..].to_vec());
            assert_eq!(curr.num_rows(), initial_rows - 1);
        }
        result.push(curr[0].clone());
        shadow_result.push(shadow[0].clone());

        *shadow = Matrix(shadow_result);

        Some(Matrix(result))
    }

    /// Perform row operations to reduce `self` to the
    /// identity matrix.  `self` must be upper triangular.
    /// All operations are performed also on `shadow` to track
    /// their cumulative effect.
    pub fn reduce_to_identity(&self, shadow: &mut Self) -> Option<Self> {
        let size = self.num_rows();
        let mut result: Vec<Vec<F>> = Vec::new();
        let mut shadow_result: Vec<Vec<F>> = Vec::new();

        for i in 0..size {
            let idx = size - i - 1;
            let row = &self.0[idx];
            let shadow_row = &shadow[idx];

            let val = row[idx];
            let inv = val.inverse()?;

            let mut normalized = scalar_vec_mul(inv, row);
            let mut shadow_normalized = scalar_vec_mul(inv, shadow_row);

            for j in 0..i {
                let idx = size - j - 1;
                let val = normalized[idx];
                let subtracted = scalar_vec_mul(val, &result[j]);
                let result_subtracted = scalar_vec_mul(val, &shadow_result[j]);

                normalized = vec_sub(&normalized, &subtracted);
                shadow_normalized =
                    vec_sub(&shadow_normalized, &result_subtracted);
            }

            result.push(normalized);
            shadow_result.push(shadow_normalized);
        }

        result.reverse();
        shadow_result.reverse();

        *shadow = Matrix(shadow_result);
        Some(Matrix(result))
    }
}

pub fn inner_product<F: PrimeField>(a: &[F], b: &[F]) -> F {
    a.iter().zip(b).fold(F::zero(), |mut acc, (v1, v2)| {
        let mut tmp = *v1;
        tmp.mul_assign(v2);
        acc.add_assign(&tmp);
        acc
    })
}

pub fn vec_add<F: PrimeField>(a: &[F], b: &[F]) -> Vec<F> {
    a.iter()
        .zip(b.iter())
        .map(|(a, b)| {
            let mut res = *a;
            res.add_assign(b);
            res
        })
        .collect::<Vec<_>>()
}

pub fn vec_sub<F: PrimeField>(a: &[F], b: &[F]) -> Vec<F> {
    a.iter()
        .zip(b.iter())
        .map(|(a, b)| {
            let mut res = *a;
            res.sub_assign(b);
            res
        })
        .collect::<Vec<_>>()
}

fn scalar_vec_mul<F: PrimeField>(scalar: F, v: &[F]) -> Vec<F> {
    v.iter()
        .map(|val| {
            let mut prod = scalar;
            prod.mul_assign(val);
            prod
        })
        .collect::<Vec<_>>()
}

pub fn kronecker_delta<F: PrimeField>(i: usize, j: usize) -> F {
    if i == j {
        F::one()
    } else {
        F::zero()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::Zero;

    type Fr = ark_bls12_381::Fr;

    #[test]
    fn test_minor() {
        let one = Fr::from(1u64);
        let two = Fr::from(2u64);
        let three = Fr::from(3u64);
        let four = Fr::from(4u64);
        let five = Fr::from(5u64);
        let six = Fr::from(6u64);
        let seven = Fr::from(7u64);
        let eight = Fr::from(8u64);
        let nine = Fr::from(9u64);

        let m: Matrix<_> = vec![
            vec![one, two, three],
            vec![four, five, six],
            vec![seven, eight, nine],
        ]
        .into();

        let cases = [
            (0, 0, Matrix(vec![vec![five, six], vec![eight, nine]])),
            (0, 1, Matrix(vec![vec![four, six], vec![seven, nine]])),
            (0, 2, Matrix(vec![vec![four, five], vec![seven, eight]])),
            (1, 0, Matrix(vec![vec![two, three], vec![eight, nine]])),
            (1, 1, Matrix(vec![vec![one, three], vec![seven, nine]])),
            (1, 2, Matrix(vec![vec![one, two], vec![seven, eight]])),
            (2, 0, Matrix(vec![vec![two, three], vec![five, six]])),
            (2, 1, Matrix(vec![vec![one, three], vec![four, six]])),
            (2, 2, Matrix(vec![vec![one, two], vec![four, five]])),
        ];
        for (i, j, expected) in &cases {
            let result = m.minor(*i, *j);

            assert_eq!(*expected, result);
        }
    }

    #[test]
    fn test_scalar_mul() {
        let zero = Fr::from(0u64);
        let one = Fr::from(1u64);
        let two = Fr::from(2u64);
        let three = Fr::from(3u64);
        let four = Fr::from(4u64);
        let six = Fr::from(6u64);

        let m = Matrix(vec![vec![zero, one], vec![two, three]]);
        let res = m.mul_by_scalar(two);

        let expected = Matrix(vec![vec![zero, two], vec![four, six]]);

        assert_eq!(expected, res);
    }

    #[test]
    fn test_vec_mul() {
        let one = Fr::from(1u64);
        let two = Fr::from(2u64);
        let three = Fr::from(3u64);
        let four = Fr::from(4u64);
        let five = Fr::from(5u64);
        let six = Fr::from(6u64);

        let a = vec![one, two, three];
        let b = vec![four, five, six];
        let res = inner_product(&a, &b);

        let expected = Fr::from(32u64);

        assert_eq!(expected, res);
    }

    #[test]
    fn test_transpose() {
        let one = Fr::from(1u64);
        let two = Fr::from(2u64);
        let three = Fr::from(3u64);
        let four = Fr::from(4u64);
        let five = Fr::from(5u64);
        let six = Fr::from(6u64);
        let seven = Fr::from(7u64);
        let eight = Fr::from(8u64);
        let nine = Fr::from(9u64);

        let m: Matrix<_> = vec![
            vec![one, two, three],
            vec![four, five, six],
            vec![seven, eight, nine],
        ]
        .into();

        let expected: Matrix<_> = vec![
            vec![one, four, seven],
            vec![two, five, eight],
            vec![three, six, nine],
        ]
        .into();

        let res = m.transpose();
        assert_eq!(expected, res);
    }

    #[test]
    fn test_inverse() {
        let zero = Fr::from(0u64);
        let one = Fr::from(1u64);
        let two = Fr::from(2u64);
        let three = Fr::from(3u64);
        let four = Fr::from(4u64);
        let five = Fr::from(5u64);
        let six = Fr::from(6u64);
        let seven = Fr::from(7u64);
        let eight = Fr::from(8u64);
        let nine = Fr::from(9u64);

        let m = Matrix(vec![
            vec![one, two, three],
            vec![four, three, six],
            vec![five, eight, seven],
        ]);

        let m1 = Matrix(vec![
            vec![one, two, three],
            vec![four, five, six],
            vec![seven, eight, nine],
        ]);

        assert!(!m1.is_invertible());
        assert!(m.is_invertible());

        let m_inv = m.invert().unwrap();

        let computed_identity = m.matmul(&m_inv).unwrap();
        assert!(computed_identity.is_identity());

        // S
        let some_vec = vec![six, five, four];

        // M^-1(S)
        let inverse_applied = m_inv.right_apply(&some_vec);

        // M(M^-1(S))
        let m_applied_after_inverse = m.right_apply(&inverse_applied);

        // S = M(M^-1(S))
        assert_eq!(
            some_vec, m_applied_after_inverse,
            "M(M^-1(V))) = V did not hold"
        );

        // panic!();
        // B
        let base_vec = vec![eight, two, five];

        // S + M(B)
        let add_after_apply = vec_add(&some_vec, &m.right_apply(&base_vec));

        // M(B + M^-1(S))
        let apply_after_add =
            m.right_apply(&vec_add(&base_vec, &inverse_applied));

        // S + M(B) = M(B + M^-1(S))
        assert_eq!(add_after_apply, apply_after_add, "breakin' the law");

        let m = Matrix(vec![vec![zero, one], vec![one, zero]]);
        let m_inv = m.invert().unwrap();
        let computed_identity = m.matmul(&m_inv).unwrap();
        assert!(computed_identity.is_identity());
        let computed_identity = m_inv.matmul(&m).unwrap();
        assert!(computed_identity.is_identity());
    }

    #[test]
    fn test_eliminate() {
        let two = Fr::from(2u64);
        let three = Fr::from(3u64);
        let four = Fr::from(4u64);
        let five = Fr::from(5u64);
        let six = Fr::from(6u64);
        let seven = Fr::from(7u64);
        let eight = Fr::from(8u64);
        let m = Matrix(vec![
            vec![two, three, four],
            vec![four, five, six],
            vec![seven, eight, eight],
        ]);

        for i in 0..m.num_rows() {
            let mut shadow = Matrix::identity(m.num_columns());
            let res = m.eliminate(i, &mut shadow);
            if i > 0 {
                assert!(res.is_none());
                continue;
            } else {
                assert!(res.is_some());
            }

            assert_eq!(
                1,
                res.unwrap()
                    .iter_rows()
                    .filter(|&row| row[i] != Fr::zero())
                    .count()
            );
        }
    }

    #[test]
    fn test_upper_triangular() {
        let zero = Fr::zero();
        let two = Fr::from(2u64);
        let three = Fr::from(3u64);
        let four = Fr::from(4u64);
        let five = Fr::from(5u64);
        let six = Fr::from(6u64);
        let seven = Fr::from(7u64);
        let eight = Fr::from(8u64);

        let m = Matrix(vec![
            vec![two, three, four],
            vec![four, five, six],
            vec![seven, eight, eight],
        ]);

        // let expected = Matrix(vec![
        //     vec![two, three, four],
        //     vec![zero, five, six],
        //     vec![zero, zero, eight],
        // ]);

        // let mut shadow = make_identity(columns(&m));
        // let _res = upper_triangular(&m, &mut shadow);
        let mut shadow = Matrix::identity(m.num_columns());
        let res = m.upper_triangular(&mut shadow).unwrap();

        // Actually assert things.
        assert!(res[0][0] != zero);
        assert!(res[0][1] != zero);
        assert!(res[0][2] != zero);
        assert!(res[1][0] == zero);
        assert!(res[1][1] != zero);
        assert!(res[1][2] != zero);
        assert!(res[2][0] == zero);
        assert!(res[2][1] == zero);
        assert!(res[2][2] != zero);
    }

    #[test]
    fn test_reduce_to_identity() {
        let two = Fr::from(2u64);
        let three = Fr::from(3u64);
        let four = Fr::from(4u64);
        let five = Fr::from(5u64);
        let six = Fr::from(6u64);
        let seven = Fr::from(7u64);
        let eight = Fr::from(8u64);

        let m = Matrix(vec![
            vec![two, three, four],
            vec![four, five, six],
            vec![seven, eight, eight],
        ]);

        let mut shadow = Matrix::identity(m.num_columns());
        let ut = m.upper_triangular(&mut shadow);

        let res = ut.and_then(|x| x.reduce_to_identity(&mut shadow)).unwrap();

        assert!(res.is_identity());

        let prod = m.matmul(&shadow).unwrap();

        assert!(prod.is_identity());
    }
}
