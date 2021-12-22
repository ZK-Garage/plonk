/// Defines a set of tests on a pairing engine / curve combination.
///
/// The set of tests is split in two. The first set between `[]` is for regular
/// tests that should not panic. The second set is for tests that should panic.
#[cfg(test)]
#[macro_export]
macro_rules! batch_test {
    ( [$($test_set:ident),*], [$($test_panic_set:ident),*] => ($engine:ty, $params:ty) ) => {
        paste::item! {
            $(
                #[test]
                #[allow(non_snake_case)]
                fn [< $test_set _on_ $engine>]() {
                    $test_set::<$engine, $params>()
                }
            )*
            $(
                #[test]
                #[should_panic]
                #[allow(non_snake_case)]
                fn [< $test_panic_set _on_ $engine>]() {
                    $test_panic_set::<$engine, $params>()
                }
            )*
        }
    }
}

#[cfg(feature = "test")]
pub use crate::constraint_system::helper::*;
