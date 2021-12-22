# PLONK

[![Build Status](https://flat.badgen.net/github/checks/zk-garage/plonk?label=ci)](https://github.com/zk-garage/plonk/actions)
[![Repository](https://flat.badgen.net/badge/icon/plonk/blue?icon=github&label=zk-garage)](https://github.com/zk-garage/plonk)

_This is a pure Rust implementation of the PLONK zk proving system_

## About

Initial implementation created by [Kev](https://github.com/kevaundray), [Carlos](https://github.com/CPerezz) and [Luke](https://github.com/LukePearson1) at Dusk Network. Redesigned by the ZK-GARAGE team to have a backend which is compatible with the [arkworks](https://github.com/arkworks-rs) suite. This allows us to leverage the multitude of curves and optimised algebra present in various arkworks repositories.

## Features

This crate includes a variety of features which will briefly be explained below:

- `parallel`: Enables `rayon` and other parallelisation primitives to be used and speed up some of the algorithms used by the crate and it's dependencies.

- `asm`: Enables inline-assembly implementations for some of the internal algorithms and primitives used by the `arkworks` dependencies of the crate.

- `trace`: Enables the Circuit debugger tooling. This is essentially the capability of using the `StandardComposer::check_circuit_satisfied` function. The function will output information about each circuit gate until one of the gates does not satisfy the equation, or there are no more gates. If there is an unsatisfied gate equation, the function will panic and return the gate number.

- `trace-print`: Goes a step further than `trace` and prints each `gate` component data, giving a clear overview of all the values which make up the circuit that we're constructing. __The recommended method is to derive the std output, and the std error, and then place them in text file which can be used to efficiently analyse the gates.__

## Documentation

There are two main types of documentation in this repository:

- **Crate documentation**. This provides info about all of the functions that the library provides, as well
  as the documentation regarding the data structures that it exports. To check this, please feel free to go to
  the [documentation page](https://docs.rs/ark-plonk/) or run `make doc` or `make doc-internal`.

- **Notes**. This is a specific subset of documentation which explains the key mathematical concepts
  of PLONK and how they work with mathematical demonstrations. To check it, run `make doc` and open the resulting docs,
  which will be located under `/target` with your browser.

## Performance

Benches taken running: `RUSTFLAGS='-C target-cpu=native' cargo bench` with an `Intel(R) Core(TM) i9-10885H`

```text
Prove  2^5 =      32 gates time:  [9.4230 ms 9.5398 ms 9.6632 ms]
Prove  2^6 =      64 gates time:  [12.896 ms 13.013 ms 13.141 ms]
Prove  2^7 =     128 gates time:  [17.836 ms 18.137 ms 18.349 ms]
Prove  2^8 =     256 gates time:  [29.553 ms 29.914 ms 30.252 ms]
Prove  2^9 =     512 gates time:  [49.285 ms 50.221 ms 51.109 ms]
Prove  2^10 =   1024 gates time:  [68.426 ms 68.704 ms 68.854 ms]
Prove  2^11 =   2048 gates time:  [127.20 ms 127.49 ms 127.70 ms]
Prove  2^12 =   4096 gates time:  [245.00 ms 245.48 ms 245.92 ms]
Prove  2^13 =   8192 gates time:  [438.23 ms 440.64 ms 442.87 ms]
Prove  2^14 =  16384 gates time:  [865.00 ms 869.66 ms 873.92 ms]
Prove  2^15 =  32768 gates time:  [1764.4 ms 1771.2 ms 1777.8 ms]
Prove  2^16 =  65536 gates time:  [3420.7 ms 3449.9 ms 3477.6 ms]
Prove  2^17 = 131072 gates time:  [6698.9 ms 6757.7 ms 6819.1 ms]
Prove  2^18 = 262144 gates time:  [13603. ms 13704. ms 13816. ms]

Verify 2^5 =      32 gates time:  [4.1958 ms 4.2881 ms 4.4535 ms]
Verify 2^6 =      64 gates time:  [4.2475 ms 4.2781 ms 4.3019 ms]
Verify 2^7 =     128 gates time:  [4.2062 ms 4.2973 ms 4.3528 ms]
Verify 2^8 =     256 gates time:  [4.1951 ms 4.2593 ms 4.3262 ms]
Verify 2^9 =     512 gates time:  [4.2094 ms 4.3023 ms 4.3846 ms]
Verify 2^10 =   1024 gates time:  [4.1785 ms 4.2228 ms 4.2789 ms]
Verify 2^11 =   2048 gates time:  [4.1177 ms 4.1379 ms 4.1461 ms]
Verify 2^12 =   4096 gates time:  [4.1342 ms 4.1467 ms 4.1632 ms]
Verify 2^13 =   8192 gates time:  [4.1329 ms 4.1770 ms 4.2117 ms]
Verify 2^14 =  16384 gates time:  [4.1817 ms 4.1870 ms 4.1921 ms]
Verify 2^15 =  32768 gates time:  [4.3285 ms 4.3390 ms 4.3591 ms]
Verify 2^16 =  65536 gates time:  [4.4808 ms 4.5020 ms 4.5267 ms]
Verify 2^17 = 131072 gates time:  [5.1339 ms 5.1572 ms 5.1806 ms]
Verify 2^18 = 262144 gates time:  [6.7577 ms 6.8124 ms 6.8925 ms]
```

## Acknowledgements

- Reference [implementation](https://github.com/AztecProtocol/barretenberg) by Aztec Protocol
- Initial [implementation](https://github.com/kobigurk/plonk/tree/kobigurk/port_to_zexe) of PLONK with arkworks backend was done years before this lib existed by Kobi Gurkan

## Licensing

This software is distributed under the terms of both the MIT license and the Apache License (Version 2.0). Please see [LICENSE-MIT](./LICENSE-MIT) and [LICENSE-APACHE](./LICENSE-APACHE) for further info.

## Contributing

- If you want to contribute to this repository/project please, check [CONTRIBUTING.md](./CONTRIBUTING.md)
- If you want to report a bug or request a new feature addition, please open an issue on this repository.
