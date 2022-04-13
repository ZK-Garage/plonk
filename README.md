# PLONK
[![CI checks](https://github.com/ZK-Garage/plonk/actions/workflows/ci.yml/badge.svg)](https://github.com/ZK-Garage/plonk/actions/workflows/ci.yml)
[![Repository](https://img.shields.io/badge/github-plonk-blueviolet?logo=github)](https://github.com/ZK-Garage/plonk)
[![Documentation](https://img.shields.io/badge/docs-plonk-blue?logo=rust)](https://docs.rs/plonk/)

_This is a pure rust implementation of various zkp components by the [ZK-Garage](https://github.com/ZK-Garage) team_

## About
This a library currently contains several modules

The `plonk-core` module is an implemention of the PLONK proving system, that leverages custom gates and lookups to significantly enhance performance and lower constraint count to optimise the generation of zero knowledge proofs. The backend of the `plonk-core` module is designed to be compatible with the [arkworks](https://github.com/arkworks-rs) suite. By leveraging the operations in arkworks we have optimised algebra and generic trait abstractions for curve types, allowing users to define their SNARK over different curves and, if applicabale, utilise embedded or pairing curves. The polynomial commitment scheme is also generic, which allows users to implement differing PCSs dependent upon particular needs.

The `plonk-hashing` module is set to contain several hashing algorithms, commencing with an optimised implementation of the Poseidon hashing algorithm generic for both plonk-style arithmetic representation and R1CS. Which will be extended but not limited to Reinforced Concrete and Blake2s.

The `plonk-book` is a module which contains a detailed overview of the working parts within a EC based ZKP system, with explanation of some of the characteristics particular to PLONK, e.g. Lagrange bases. There is a also a chapter on the construction of the PLONK algorithms, as well as an explanation of the implementations features and details specific to this repository.
### Compile the Plonk book
First, you need to install mdbook command line tool used to create books with Markdown.

```cargo install mdbook```



You should also install Katex preprocessor which renders Latex equations into HTML at build time

```cargo install --git "https://github.com/lzanini/mdbook-katex"```


Then, you build the book as follows:

```mdbook build```

Last but not least, you can read the book by doing this command

```mdbook serve --open```

This will display the book in your default web browser after building it.



### Join the effort

Please, if you're interested in collaborating, contributing or just discussing, you can join our Discord here: <https://discord.gg/XWJdhVf37F>

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
  which will be located under `/target/doc/plonk/index.html` with your browser.

  **Examples**. Examples can be found in the `examples` folder. Run them, e.g., via `cargo run --example simple_circuit`.

## Performance

Benches taken running: `RUSTFLAGS='-C target-cpu=native' cargo bench` with an `AMD Ryzen 7 3700X`
These benches use the `Bls12-381` curve.

Using `KZG10` commitments:
```
Compile:
2^5         [    17.632 ms      17.669 ms      17.696 ms]
2^6         [    22.666 ms      22.702 ms      22.747 ms]
2^7         [    29.618 ms      29.653 ms      29.719 ms]
2^8         [    47.467 ms      47.556 ms      47.609 ms]
2^9         [    65.458 ms      65.786 ms      66.174 ms]
2^10        [    97.172 ms      97.514 ms      97.897 ms]
2^11        [   167.89  ms     168.17  ms     168.41  ms]
2^12        [   314.51  ms     314.65  ms     314.78  ms]
2^13        [   526.59  ms     527.63  ms     529.18  ms]
2^14        [ 1.0238     s   1.0253     s   1.0272     s]
2^15        [ 2.0029     s   2.0088     s   2.0143     s]
2^16        [ 3.7727     s   3.7846     s   3.7955     s]
2^17        [ 6.7340     s   6.7523     s   6.7700     s]
2^18        [13.584      s  13.613      s  13.640      s]

Prove:
2^5         [    16.172 ms      16.208 ms      16.264 ms]
2^6         [    21.676 ms      21.712 ms      21.748 ms]
2^7         [    29.493 ms      29.545 ms      29.613 ms]
2^8         [    48.970 ms      49.039 ms      49.104 ms]
2^9         [    72.251 ms      72.533 ms      72.703 ms]
2^10        [   128.89  ms     130.71  ms     132.10  ms]
2^11        [   242.91  ms     247.74  ms     252.29  ms]
2^12        [   455.79  ms     459.45  ms     462.85  ms]
2^13        [   776.94  ms     781.89  ms     787.94  ms]
2^14        [ 1.4752     s   1.4824     s   1.4893     s]
2^15        [ 2.8589     s   2.8682     s   2.8787     s]
2^16        [ 5.4610     s   5.4766     s   5.4927     s]
2^17        [10.078      s  10.118      s  10.159      s]
2^18        [20.151      s  20.184      s  20.216      s]

Verify:
2^5         [    5.5250 ms      5.5560 ms      5.5983 ms]
2^6         [    5.4933 ms      5.5461 ms      5.5910 ms]
2^7         [    5.5678 ms      5.6002 ms      5.6247 ms]
2^8         [    5.5391 ms      5.5756 ms      5.6027 ms]
2^9         [    5.5421 ms      5.5648 ms      5.5929 ms]
2^10        [    5.5423 ms      5.5825 ms      5.6240 ms]
2^11        [    5.5269 ms      5.5576 ms      5.6027 ms]
2^12        [    5.5624 ms      5.6081 ms      5.6623 ms]
2^13        [    5.6288 ms      5.6656 ms      5.6914 ms]
2^14        [    5.6068 ms      5.6186 ms      5.6292 ms]
2^15        [    5.5930 ms      5.6241 ms      5.6543 ms]
2^16        [    6.0845 ms      6.1324 ms      6.1745 ms]
2^17        [    6.5760 ms      6.5896 ms      6.6030 ms]
2^18        [    8.1152 ms      8.1481 ms      8.1710 ms]
```
Using `IPA` commitments:
```
Compile:
2^5         [    16.768 ms      16.818 ms      16.857 ms]
2^6         [    21.958 ms      21.977 ms      21.993 ms]
2^7         [    28.847 ms      28.869 ms      28.903 ms]
2^8         [    47.626 ms      47.660 ms      47.693 ms]
2^9         [    67.319 ms      67.485 ms      67.674 ms]
2^10        [    98.526 ms      98.891 ms      99.072 ms]
2^11        [   171.84  ms     172.06  ms     172.25  ms]
2^12        [   322.42  ms     322.55  ms     322.69  ms]
2^13        [   533.50  ms     533.95  ms     534.53  ms]
2^14        [ 1.0333     s   1.0342     s   1.0351     s]
2^15        [ 2.0156     s   2.0240     s   2.0308     s]
2^16        [ 3.8668     s   3.8769     s   3.8871     s]
2^17        [ 6.8066     s   6.8259     s   6.8506     s]
2^18        [13.757      s  13.773      s  13.788      s]

Prove:
2^5         [    32.205 ms      32.802 ms      33.418 ms]
2^6         [    39.419 ms      39.479 ms      39.550 ms]
2^7         [    53.665 ms      53.767 ms      53.876 ms]
2^8         [    83.829 ms      84.005 ms      84.171 ms]
2^9         [   127.58  ms     127.85  ms     128.11  ms]
2^10        [   207.01  ms     208.50  ms     210.09  ms]
2^11        [   397.91  ms     400.53  ms     403.63  ms]
2^12        [   719.49  ms     725.85  ms     732.68  ms]
2^13        [ 1.2864     s   1.2912     s   1.2953     s]
2^14        [ 2.4494     s   2.4552     s   2.4620     s]
2^15        [ 4.7411     s   4.7617     s   4.7826     s]
2^16        [ 9.1925     s   9.2148     s   9.2360     s]
2^17        [17.499      s  17.584      s  17.660      s]
2^18        [35.019      s  35.084      s  35.138      s]

Verify:
2^5         [    7.9861 ms      8.0159 ms      8.0433 ms]
2^6         [    8.9787 ms      9.0031 ms      9.0272 ms]
2^7         [   10.648  ms     10.675  ms     10.714  ms]
2^8         [   13.466  ms     13.526  ms     13.596  ms]
2^9         [   17.140  ms     17.188  ms     17.267  ms]
2^10        [   25.379  ms     25.574  ms     25.785  ms]
2^11        [   34.424  ms     37.413  ms     38.720  ms]
2^12        [   39.254  ms     39.429  ms     39.595  ms]
2^13        [   69.872  ms     70.392  ms     70.790  ms]
2^14        [  130.16   ms    130.93   ms    131.90   ms]
2^15        [  243.71   ms    246.59   ms    249.40   ms]
2^16        [  409.56   ms    415.00   ms    419.81   ms]
2^17        [  777.07   ms    789.39   ms    801.28   ms]
2^18        [1.4931      s  1.4999      s  1.5065      s]
```

## Acknowledgements

- Reference [implementation](https://github.com/AztecProtocol/barretenberg) by Aztec Protocol
- Initial [implementation](https://github.com/kobigurk/plonk/tree/kobigurk/port_to_zexe) of PLONK with arkworks backend was done years before this lib existed by Kobi Gurkan
- Initial rust [implementation](https://github.com/dusk-network/plonk) of PLONK by Dusk Network, originally forked from dusk-network/plonk

## Licensing

This software is distributed under the terms of Mozilla Public License Version 2.0 (MPL-2.0). Please see [LICENSE](https://github.com/ZK-Garage/plonk/blob/master/LICENSE) for further info.

## Contributing
- If you want to contribute to this repository/project please, check [CONTRIBUTING.md](./CONTRIBUTING.md)
- If you want to report a bug or request a new feature addition, please open an issue on this repository.
