// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright (c) ZK-GARAGE. All rights reserved.

//! Test Suite

/// Defines a set of tests on a pairing engine / curve combination.
///
/// The set of tests is split in two. The first set between `[]` is for regular
/// tests that should not panic. The second set is for tests that should panic.

#[macro_export]
macro_rules! batch_test_field {
    ( [$($test_set:ident),*], [$($test_panic_set:ident),*] => ($engine:ty) ) => {
        paste::item! {
            $(
                #[test]
                #[allow(non_snake_case)]
                fn [< $test_set _on_ $engine>]() {
                    $test_set::<<$engine as ark_ec::PairingEngine>::Fr>()
                }
            )*
            $(
                #[test]
                #[should_panic]
                #[allow(non_snake_case)]
                fn [< $test_panic_set _on_ $engine>]() {
                    $test_panic_set::<<$engine as ark_ec::PairingEngine>::Fr>()
                }
            )*
        }
    }
}

#[macro_export]
macro_rules! batch_test_embedded {

    ([$($test_set:ident),*], [$($test_panic_set:ident),*] => [$param:ty, $($params:ty),+]) => {
        batch_test_embedded! { [$($test_set),*], [$($test_panic_set),*] => [$param] }
        batch_test_embedded! { [$($test_set),*], [$($test_panic_set),*] => [$($params),+] }
    };

     ([$($test_set:ident),*], [$($test_panic_set:ident),*] => [$params:ty])  => {
        paste::item! {
            $(
                #[test]
                #[allow(non_snake_case)]
                fn [< $test_set _on_ $params >]() {

                    $test_set::<$params, _, _>()
                }
            )*
            $(
                #[test]
                #[should_panic]
                #[allow(non_snake_case)]
                fn [< $test_panic_set _on_ $params >]() {
                    $test_panic_set::<$params, _, _>()
                }
            )*
        }
    };
}

#[macro_export]
macro_rules! batch_test {

    ([$($test_set:ident),*], [$($test_panic_set:ident),*] => [$param:ty, $($params:ty),+]) => {
        batch_test! { [$($test_set),*], [$($test_panic_set),*] => [$param] }
        batch_test! { [$($test_set),*], [$($test_panic_set),*] => [$($params),+] }
    };

     ([$($test_set:ident),*], [$($test_panic_set:ident),*] => [$params:ty])  => {
        paste::item! {
            $(
                #[test]
                #[allow(non_snake_case)]
                fn [< $test_set _on_ $params >]() {

                    $test_set::<$params>()
                }
            )*
            $(
                #[test]
                #[should_panic]
                #[allow(non_snake_case)]
                fn [< $test_panic_set _on_ $params >]() {
                    $test_panic_set::<$params>()
                }
            )*
        }
    };
}

#[macro_export]
macro_rules! batch_field_test {
    ( [$($test_set:ident),*], [$($test_panic_set:ident),*] => $field:ty ) => {
        paste::item! {
            $(
                #[test]
                #[allow(non_snake_case)]
                fn [< $test_set _on_ $field>]() {
                    $test_set::<$field>()
                }
            )*
            $(
                #[test]
                #[should_panic]
                #[allow(non_snake_case)]
                fn [< $test_panic_set _on_ $field>]() {
                    $test_panic_set::<$field>()
                }
            )*
        }
    }
}
