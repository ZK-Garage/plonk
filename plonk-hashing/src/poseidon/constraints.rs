//! Library independent specification for field trait.
//! `COM` is `()` when the field is in native, and is constraint synthesizer
//! when the field is a variable.
//!
//! This file is adapted from manta-rs, but is tweaked to better support PLONK proving system.
//! Some essential features (like allocating unknown variable) are missing.
//!
//! The initial goal is to simultaneous support PLONK, R1CS-std, and native.

/// Compiler constant.
pub trait Constant<COM>
where COM: ?Sized
{
    type Type;
    fn new_constant(c: &mut COM, value: &Self::Type) -> Self;
}
/// Constant Type Alias
pub type Const<C, COM> = <C as Constant<COM>>::Type;

/// Compiler variable.
pub trait Variable<COM>
where COM: ?Sized {
    type Type;
    // we can only make it unknown (the witness), because PLONK does not support arbitrary public input.
    fn new_variable(c: &mut COM, value: &Self::Type) -> Self;
}
/// Variable Type Alias
pub type Var<V, COM> = <V as Variable<COM>>::Type;

/// Value-source auto-trait
pub trait ValueSource<COM>
where COM: ?Sized
{
    /// Allocates `self` as a constant in `compiler`.
    fn as_constant<C>(&self, compiler: &mut COM) -> C
        where
            C: Constant<COM, Type = Self>,
    {
        C::new_constant(compiler, self)
    }

    /// Allocates `self` as a known value in `compiler`.
    fn as_variable<M, V>(&self, compiler: &mut COM) -> V
        where
            V: Variable<COM, Type = Self>,
    {
        V::new_variable(compiler, self)
    }
}

impl<COM, T> ValueSource<COM> for T where T: ?Sized {}

/// Allocator Auto-Trait
pub trait Allocator {
    /// Allocates a constant with the given `value` into `self`.
    #[inline]
    fn allocate_constant<C>(&mut self, value: &C::Type) -> C
        where
            C: Constant<Self>,
    {
        C::new_constant(self, value)
    }

    /// Allocates a known variable with the given `value` into `self`.
    #[inline]
    fn allocate_variable<M, V>(&mut self, value: &V::Type) -> V
        where
            V: Variable<Self>,
    {
        V::new_variable(self, value)
    }
}

impl<COM> Allocator for COM where COM: ?Sized {}

/// Addition trait for variable
pub trait COMAdd<COM>: Variable<COM>
    where
        COM: ?Sized,
{
    type Constant: Constant<COM>;
    /// Adds `self` and `rhs` inside of `compiler`.
    fn com_add(&self, rhs: &Self, compiler: &mut COM) -> Self;
}

/// Subtraction
pub trait COMSub<COM>: Sized
    where
        COM: ?Sized,
{
    type Output;
    /// Subtracts `rhs` from `self` inside of `compiler`.
    fn com_sub(&self, rhs: &Self, compiler: &mut COM) -> Self;
}

/// Constraint System trait for compiler
pub trait ConstraintSystem {
    /// Boolean Variable type
    type Bool;

    /// Asserts that `b == 1`.
    fn assert(&mut self, b: Self::Bool);

    /// Asserts that all the booleans in `iter` are equal to `1`.
    #[inline]
    fn assert_all<I>(&mut self, iter: I)
        where
            I: IntoIterator<Item = Self::Bool>,
    {
        iter.into_iter().for_each(move |b| self.assert(b));
    }

    /// Generates a boolean that represents the fact that `lhs` and `rhs` may be equal.
    #[inline]
    fn eq<V>(&mut self, lhs: &V, rhs: &V) -> Self::Bool
        where
            V: Equal<Self>,
    {
        V::eq(lhs, rhs, self)
    }

    /// Asserts that `lhs` and `rhs` are equal.
    #[inline]
    fn assert_eq<V>(&mut self, lhs: &V, rhs: &V)
        where
            V: Equal<Self>,
    {
        V::assert_eq(lhs, rhs, self);
    }
}

/// Equality Trait
pub trait Equal<COM>
    where
        COM: ConstraintSystem + ?Sized,
{
    /// Generates a boolean that represents the fact that `lhs` and `rhs` may be equal.
    fn eq(lhs: &Self, rhs: &Self, compiler: &mut COM) -> COM::Bool;

    /// Asserts that `lhs` and `rhs` are equal.
    #[inline]
    fn assert_eq(lhs: &Self, rhs: &Self, compiler: &mut COM) {
        let boolean = Self::eq(lhs, rhs, compiler);
        compiler.assert(boolean);
    }

    /// Asserts that all the elements in `iter` are equal to some `base` element.
    #[inline]
    fn assert_all_eq_to_base<'t, I>(base: &'t Self, iter: I, compiler: &mut COM)
        where
            I: IntoIterator<Item = &'t Self>,
    {
        for item in iter {
            Self::assert_eq(base, item, compiler);
        }
    }

    /// Asserts that all the elements in `iter` are equal.
    #[inline]
    fn assert_all_eq<'t, I>(iter: I, compiler: &mut COM)
        where
            Self: 't,
            I: IntoIterator<Item = &'t Self>,
    {
        let mut iter = iter.into_iter();
        if let Some(base) = iter.next() {
            Self::assert_all_eq_to_base(base, iter, compiler);
        }
    }

}




// /// Basic Arithmetic Operations. Both constants and variables should implement this trait.
// pub trait COMArith<COM = ()>: Sized + Clone + Debug {
//     // I added `com_` prefix here to avoid conflict with num_traits. Any suggestion
//     // is welcome!
//     // TODO: do we want to split `add`, `neg`, `mul` to different traits?
//     /// additive identity
//     fn com_zero(c: &mut COM) -> Self;
//     fn zeros<const SIZE: usize>(c: &mut COM) -> [Self; SIZE];
//     /// add two field elements
//     fn com_add(&self, c: &mut COM, b: &Self) -> Self;
//     /// the additive inverse of a field element
//     fn com_neg(&self, c: &mut COM) -> Self;
//     /// multiply two field elements
//     fn com_mul(&self, c: &mut COM, other: &Self) -> Self;
//     fn com_square(&self, c: &mut COM) -> Self {
//         self.com_mul(c, self)
//     }
//     fn com_add_assign(&mut self, c: &mut COM, other: &Self) {
//         *self = self.com_add(c, other);
//     }
//     fn com_mul_assign(&mut self, c: &mut COM, other: &Self) {
//         *self = self.com_mul(c, other);
//     }
// }
//
// pub trait COMFromInt<COM = ()>: COMArith<COM> {
//     /// multiplicative identity
//     fn com_one(c: &mut COM) -> Self;
//     fn com_from_const_int(c: &mut COM, v: u64) -> Self;
// }
//
// pub trait COMPower<COM = ()>: COMArith<COM> {
//     type Scalar;
//     fn com_pow(&self, c: &mut COM, exp: &Self::Scalar) -> Self;
// }
//
// pub trait COMArithExt<COM = ()>: COMArith<COM> + Sized {
//     fn __make_arith_gate(c: &mut COM, config: ArithExtBuilder<Self, COM>) -> Self;
//     /// `(w_l * w_r) * q_m + a * q_l + b * q_r + w_4 * q_4 + q_c + PI + q_o * c = 0`
//     /// where output is `c`
//     fn com_arith(c: &mut COM) -> ArithExtBuilder<Self, COM> {
//         ArithExtBuilder::new(c)
//     }
// }
//
// pub type NativeField<S, COM> = <S as COMArith<COM>>::Native;
//
// /// `(w_l * w_r) * q_m + a * q_l + b * q_r + w_4 * q_4 + q_c + PI + q_o * c = 0`
// /// where output is `c`
// pub struct ArithExtBuilder<F: COMArithExt<COM>, COM = ()> {
//     w_l: F,
//     w_r: F,
//     q_m: F::Native,
//     q_l: F::Native,
//     q_r: F::Native,
//     q_c: F::Native,
//     q_o: F::Native,
//     q_4_w_4: Option<(F::Native, F)>,
//     pi: Option<F::PublicInput>,
//     _compiler: PhantomData<COM>,
// }
//
// impl<F: COMArithExt<COM>, COM> ArithExtBuilder<F, COM> {
//     pub(crate) fn new(c: &mut COM) -> Self {
//         Self {
//             w_l: F::com_zero(c),
//             w_r: F::com_zero(c),
//             q_m: F::Native::com_zero(&mut ()),
//             q_l: F::Native::com_zero(&mut ()),
//             q_r: F::Native::com_zero(&mut ()),
//             q_c: F::Native::com_zero(&mut ()),
//             q_o: F::Native::com_one(&mut ()).com_neg(&mut ()),
//             q_4_w_4: None,
//             pi: None,
//             _compiler: PhantomData,
//         }
//     }
//
//     pub fn w_l(mut self, w_l: F) -> Self {
//         self.w_l = w_l;
//         self
//     }
//
//     pub fn w_r(mut self, w_r: F) -> Self {
//         self.w_r = w_r;
//         self
//     }
//
//     pub fn witness(mut self, w_l: F, w_r: F) -> Self {
//         self.w_l = w_l;
//         self.w_r = w_r;
//         self
//     }
//
//     pub fn q_m(mut self, q_m: F::Native) -> Self {
//         self.q_m = q_m;
//         self
//     }
//
//     pub fn mul(mut self) -> Self {
//         self.q_m(F::Native::com_one(&mut ()))
//     }
//
//     pub fn q_l(mut self, q_l: F::Native) -> Self {
//         self.q_l = q_l;
//         self
//     }
//
//     pub fn q_r(mut self, q_r: F::Native) -> Self {
//         self.q_r = q_r;
//         self
//     }
//
//     pub fn q_c(mut self, q_c: F::Native) -> Self {
//         self.q_c = q_c;
//         self
//     }
//
//     pub fn q4w4(mut self, q_4_w_4: (F::Native, F)) -> Self {
//         self.q_4_w_4 = Some(q_4_w_4);
//         self
//     }
//
//     pub fn pi(mut self, pi: F::PublicInput) -> Self {
//         self.pi = Some(pi);
//         self
//     }
//
//     pub fn q_o(mut self, q_o: F::Native) -> Self {
//         self.q_o = q_o;
//         self
//     }
//
//     pub fn build(self, c: &mut COM) -> F {
//         F::__make_arith_gate(c, self)
//     }
// }
//
// impl<F: PrimeField> COMArith<()> for F {
//     type Native = F;
//     type PublicInput = F;
//     fn com_zero(_c: &mut ()) -> Self {
//         F::zero()
//     }
//
//     fn com_alloc(_c: &mut (), v: Self::Native) -> Self {
//         v
//     }
//
//     fn zeros<const SIZE: usize>(_c: &mut ()) -> [Self; SIZE] {
//         [F::zero(); SIZE]
//     }
//
//     fn com_add(&self, _c: &mut (), b: &Self) -> Self {
//         *self + *b
//     }
//
//     fn com_addi(&self, _c: &mut (), b: &Self::Native) -> Self {
//         *self + *b
//     }
//
//     fn com_neg(&self, _c: &mut ()) -> Self {
//         -*self
//     }
//
//     fn com_mul(&self, _c: &mut (), other: &Self) -> Self {
//         *self * *other
//     }
//
//     fn com_muli(&self, _c: &mut (), other: &Self::Native) -> Self {
//         *self * *other
//     }
// }
//
// impl<F: PrimeField> COMFromInt<()> for F {
//     fn com_one(_c: &mut ()) -> Self {
//         F::one()
//     }
//
//     fn com_from_const_int(c: &mut (), v: u64) -> Self {
//         F::from(v)
//     }
// }
//
// impl<F: PrimeField> COMArithExt<()> for F {
//     fn __make_arith_gate(_c: &mut (), config: ArithExtBuilder<Self, ()>) -> Self {
//         let mut result = F::zero();
//         result += (config.w_l * config.w_r) * config.q_m;
//         result += config.q_l * config.w_l;
//         result += config.q_r * config.w_r;
//         result += config.q_4_w_4.map_or(F::zero(), |(q_4, w_4)| q_4 * w_4);
//         result += config.q_c;
//         result += config.pi.unwrap_or(F::zero());
//
//         // now result = - q_o * c, we want c = (-result) / q_o
//         let q_o_inv = F::inverse(&config.q_o).unwrap();
//         (-result) * q_o_inv
//     }
// }
//
// impl<E, P> COMArith<StandardComposer<E, P>> for Variable
// where
//     E: PairingEngine,
//     P: TEModelParameters<BaseField = E::Fr>,
// {
//     type Native = E::Fr;
//     type PublicInput = E::Fr;
//
//     fn com_zero(c: &mut StandardComposer<E, P>) -> Self {
//         c.zero_var()
//     }
//
//     fn com_alloc(c: &mut StandardComposer<E, P>, v: Self::Native) -> Self {
//         c.add_input(v)
//     }
//
//     fn zeros<const SIZE: usize>(c: &mut StandardComposer<E, P>) -> [Self; SIZE] {
//         [c.zero_var(); SIZE]
//     }
//
//     fn com_add(&self, c: &mut StandardComposer<E, P>, b: &Self) -> Self {
//         c.arithmetic_gate(|g| g.witness(*self, *b, None).add(E::Fr::one(), E::Fr::one()))
//     }
//
//     fn com_addi(&self, c: &mut StandardComposer<E, P>, b: &Self::Native) -> Self {
//         let zero = c.zero_var();
//         c.arithmetic_gate(|g| {
//             g.witness(*self, zero, None)
//                 .add(E::Fr::one(), E::Fr::zero())
//                 .constant(*b)
//         })
//     }
//
//     fn com_neg(&self, c: &mut StandardComposer<E, P>) -> Self {
//         let zero = c.zero_var();
//         c.arithmetic_gate(|g| g.witness(*self, zero, None).out(E::Fr::one()))
//     }
//
//     fn com_mul(&self, c: &mut StandardComposer<E, P>, other: &Self) -> Self {
//         c.arithmetic_gate(|g| g.witness(*self, *other, None).mul(E::Fr::one()))
//     }
//
//     fn com_muli(&self, c: &mut StandardComposer<E, P>, other: &Self::Native) -> Self {
//         let zero = c.zero_var();
//         c.arithmetic_gate(|g| g.witness(*self, zero, None).add(*other, E::Fr::zero()))
//     }
// }
//
// impl<E, P> COMArithExt<StandardComposer<E, P>> for Variable
// where
//     E: PairingEngine,
//     P: TEModelParameters<BaseField = E::Fr>,
// {
//     fn __make_arith_gate(
//         c: &mut StandardComposer<E, P>,
//         config: ArithExtBuilder<Self, StandardComposer<E, P>>,
//     ) -> Self {
//         c.arithmetic_gate(|g| {
//             g.witness(config.w_l, config.w_r, None).mul(config.q_m);
//             g.add(config.q_l, config.q_r);
//             if let Some((q_4, w_4)) = config.q_4_w_4 {
//                 g.fan_in_3(q_4, w_4);
//             };
//             g.constant(config.q_c);
//             g.out(config.q_o);
//             if let Some(pi) = config.pi {
//                 g.constant(pi);
//             }
//             g
//         })
//     }
// }
//
// #[cfg(test)]
// mod tests {
//     use super::*;
//     use ark_bls12_381::Fr;
//     use ark_std::{test_rng, UniformRand};
//
//     #[test]
//     fn sanity_check_on_native() {
//         // calculate 3xy + 2x + y + 1
//         let mut rng = test_rng();
//         let x = Fr::rand(&mut rng);
//         let y = Fr::rand(&mut rng);
//         let expected = (Fr::from(3u64) * x * y) + (Fr::from(2u64) * x) + y + Fr::one();
//         let actual = Fr::com_arith(&mut ())
//             .w_l(x)
//             .w_r(y)
//             .q_m(3u64.into())
//             .q_l(2u64.into())
//             .q_r(Fr::one())
//             .q_c(Fr::one())
//             .build(&mut ());
//
//         assert_eq!(expected, actual);
//     }
// }
