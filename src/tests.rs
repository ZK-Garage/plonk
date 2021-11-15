#[cfg(test)]
mod test {
    /// Run a set of tests on a pairing engine / curve combination
    #[macro_export]
    macro_rules! batch_test {
        ( [$($test_set:ident),+] => ($engine:ty, $curve:ty, $params:ty) ) => {
            paste::item! {
                $(
                #[test]
                #[allow(non_snake_case)]
                fn [< $test_set _on_ $engine>]() {
                    $test_set::<$engine, $curve, $params>()
                }
                    )+
            }
        }

    }
}
