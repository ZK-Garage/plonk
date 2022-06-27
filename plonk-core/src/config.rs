use crate::commitment::HomomorphicCommitment;

trait Configuration {
    /// Total Number of Wires Accessible to the Compiler
    const WIRE_COUNT: usize;

    /// Scalar Field Type
    ///
    /// This type represents the underlying coefficient field for polynomial
    /// constraints.
    type Field: ark_ff::PrimeField;

    /// Polynomial Commitment Scheme Type
    type PolynomialCommitment: HomomorphicCommitment<Self::Field>;

    /// Polynomial Commitment Opening Type
    type Opening;

    //TODO Complete with all necessary fields
}
