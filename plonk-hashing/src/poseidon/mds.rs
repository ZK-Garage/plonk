// adapted from https://github.com/filecoin-project/neptune/blob/master/src/mds.rs
use crate::poseidon::matrix::Matrix;
use ark_ff::vec::*;
use ark_ff::PrimeField;

#[derive(Clone, Debug, PartialEq)]
pub struct MdsMatrices<F: PrimeField> {
    pub m: Matrix<F>,
    pub m_inv: Matrix<F>,
    pub m_hat: Matrix<F>,
    pub m_hat_inv: Matrix<F>,
    pub m_prime: Matrix<F>,
    pub m_double_prime: Matrix<F>,
}

impl<F: PrimeField> MdsMatrices<F> {
    /// Derive MDS matrix of size `dim*dim` and relevant things
    pub fn new(dim: usize) -> Self {
        let m = Self::generate_mds(dim);
        Self::derive_mds_matrices(m)
    }

    /// Given an MDS matrix `m`, compute all its associated matrices.
    pub(crate) fn derive_mds_matrices(m: Matrix<F>) -> Self {
        let m_inv = m.invert().expect("Derived MDS matrix is not invertible");
        let m_hat = m.minor(0, 0);
        let m_hat_inv =
            m_hat.invert().expect("Derived MDS matrix is not correct");
        let m_prime = Self::make_prime(&m);
        let m_double_prime = Self::make_double_prime(&m, &m_hat_inv);
        MdsMatrices {
            m,
            m_inv,
            m_hat,
            m_hat_inv,
            m_prime,
            m_double_prime,
        }
    }

    fn generate_mds(t: usize) -> Matrix<F> {
        let xs: Vec<F> = (0..t as u64).map(F::from).collect();
        let ys: Vec<F> = (t as u64..2 * t as u64).map(F::from).collect();

        let matrix = xs
            .iter()
            .map(|xs_item| {
                ys.iter()
                    .map(|ys_item| {
                        // Generate the entry at (i,j)
                        let mut tmp = *xs_item;
                        tmp.add_assign(ys_item);
                        tmp.inverse().unwrap()
                    })
                    .collect()
            })
            .collect::<Matrix<F>>();

        debug_assert!(matrix.is_invertible());
        debug_assert_eq!(matrix, matrix.transpose());
        matrix
    }

    /// Returns a matrix associated to `m` in the optimization of
    /// MDS matrices.
    fn make_prime(m: &Matrix<F>) -> Matrix<F> {
        m.iter_rows()
            .enumerate()
            .map(|(i, row)| match i {
                0 => {
                    let mut new_row = vec![F::zero(); row.len()];
                    new_row[0] = F::one();
                    new_row
                }
                _ => {
                    let mut new_row = vec![F::zero(); row.len()];
                    new_row[1..].copy_from_slice(&row[1..]);
                    new_row
                }
            })
            .collect()
    }

    /// Returns a matrix associated to `m` in the optimization of
    /// MDS matrices.
    fn make_double_prime(m: &Matrix<F>, m_hat_inv: &Matrix<F>) -> Matrix<F> {
        let (v, w) = Self::make_v_w(m);
        let w_hat = m_hat_inv.right_apply(&w);

        m.iter_rows()
            .enumerate()
            .map(|(i, row)| match i {
                0 => {
                    let mut new_row = Vec::with_capacity(row.len());
                    new_row.push(row[0]);
                    new_row.extend(&v);
                    new_row
                }
                _ => {
                    let mut new_row = vec![F::zero(); row.len()];
                    new_row[0] = w_hat[i - 1];
                    new_row[i] = F::one();
                    new_row
                }
            })
            .collect()
    }

    /// Returns two vectors associated to `m` in the optimization of
    /// MDS matrices.
    fn make_v_w(m: &Matrix<F>) -> (Vec<F>, Vec<F>) {
        let v = m[0][1..].to_vec();
        let w = m.iter_rows().skip(1).map(|column| column[0]).collect();
        (v, w)
    }
}

/// A `SparseMatrix` is specifically one of the form of M''.
/// This means its first row and column are each dense, and the interior matrix
/// (minor to the element in both the row and column) is the identity.
/// We will pluralize this compact structure `sparse_matrixes` to distinguish
/// from `sparse_matrices` from which they are created.
#[derive(Debug, Clone, PartialEq)]
pub struct SparseMatrix<F: PrimeField> {
    /// `w_hat` is the first column of the M'' matrix. It will be directly
    /// multiplied (scalar product) with a row of state elements.
    pub w_hat: Vec<F>,
    /// `v_rest` contains all but the first (already included in `w_hat`).
    pub v_rest: Vec<F>,
}

impl<F: PrimeField> SparseMatrix<F> {
    pub fn new(m_double_prime: &Matrix<F>) -> Self {
        assert!(m_double_prime.is_sparse());

        let w_hat = m_double_prime.iter_rows().map(|r| r[0]).collect();
        let v_rest = m_double_prime[0][1..].to_vec();
        Self { w_hat, v_rest }
    }

    pub fn size(&self) -> usize {
        self.w_hat.len()
    }

    pub fn to_matrix(&self) -> Matrix<F> {
        let mut m = Matrix::identity(self.size());
        for (j, elt) in self.w_hat.iter().enumerate() {
            m[j][0] = *elt;
        }
        for (i, elt) in self.v_rest.iter().enumerate() {
            m[0][i + 1] = *elt;
        }
        m
    }
}

// TODO: naming is from https://github.com/filecoin-project/neptune/blob/master/src/mds.rs
// TODO: which is little difficult to follow... We need to change it

pub fn factor_to_sparse_matrixes<F: PrimeField>(
    base_matrix: Matrix<F>,
    n: usize,
) -> (Matrix<F>, Vec<SparseMatrix<F>>) {
    let (pre_sparse, mut sparse_matrices) =
        (0..n).fold((base_matrix.clone(), Vec::new()), |(curr, mut acc), _| {
            let derived = MdsMatrices::derive_mds_matrices(curr);
            acc.push(derived.m_double_prime);
            let new = base_matrix.matmul(&derived.m_prime).unwrap();
            (new, acc)
        });
    sparse_matrices.reverse();
    let sparse_matrixes = sparse_matrices
        .iter()
        .map(|m| SparseMatrix::<F>::new(m))
        .collect::<Vec<_>>();

    (pre_sparse, sparse_matrixes)
}

#[cfg(test)]
mod tests {
    use crate::poseidon::mds::MdsMatrices;
    use ark_bls12_381::Fr;
    use ark_std::{test_rng, UniformRand};

    #[test]
    fn test_mds_matrices_creation() {
        for i in 2..5 {
            test_mds_matrices_creation_aux(i);
        }
    }

    fn test_mds_matrices_creation_aux(width: usize) {
        let MdsMatrices {
            m,
            m_inv,
            m_hat,
            m_hat_inv: _,
            m_prime,
            m_double_prime,
        } = MdsMatrices::<Fr>::new(width);

        for i in 0..m_hat.num_rows() {
            for j in 0..m_hat.num_columns() {
                assert_eq!(
                    m[i + 1][j + 1],
                    m_hat[i][j],
                    "MDS minor has wrong value."
                );
            }
        }

        // M^-1 x M = I
        assert!(m_inv.matmul(&m).unwrap().is_identity());

        // M' x M'' = M
        assert_eq!(m, m_prime.matmul(&m_double_prime).unwrap());
    }

    #[test]
    fn test_swapping() {
        test_swapping_aux(3)
    }

    fn test_swapping_aux(width: usize) {
        let mut rng = test_rng();
        let mds = MdsMatrices::<Fr>::new(width);

        let base = (0..width).map(|_| Fr::rand(&mut rng)).collect::<Vec<_>>();
        let x = {
            let mut x = base.clone();
            x[0] = Fr::rand(&mut rng);
            x
        };
        let y = {
            let mut y = base;
            y[0] = Fr::rand(&mut rng);
            y
        };

        let qx = mds.m_prime.right_apply(&x);
        let qy = mds.m_prime.right_apply(&y);
        assert_eq!(qx[0], x[0]);
        assert_eq!(qy[0], y[0]);
        assert_eq!(qx[1..], qy[1..]);

        let mx = mds.m.left_apply(&x);
        let m1_m2_x =
            mds.m_prime.left_apply(&mds.m_double_prime.left_apply(&x));
        assert_eq!(mx, m1_m2_x);

        let xm = mds.m.right_apply(&x);
        let x_m1_m2 =
            mds.m_double_prime.right_apply(&mds.m_prime.right_apply(&x));
        assert_eq!(xm, x_m1_m2);
    }
}
