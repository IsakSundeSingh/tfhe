# TFHE

[![Continuous Integration](https://github.com/IsakSundeSingh/tfhe/workflows/Continuous%20Integration/badge.svg)](https://github.com/IsakSundeSingh/tfhe/actions)

TFHE is a Rust port of the C++ library [TFHE](https://github.com/tfhe/tfhe).

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
tfhe = { git = "https://github.com/IsakSundeSingh/tfhe" }
```

## Documentation

To generate and open the user-documentation for the crate in your browser, clone the repository and run the following:

```bash
cargo doc --all-features --open
```

If you wish to contribute to the code and want internal documentation for implementation details, run the following:

```bash
cargo doc --all-features --document-private-items --open
```

## Crate Features

TFHE is built with these features enabled by default:

- `fft` enables FFT for polynomial multiplication, improving speeds. No reason to disable this unless you wish to remove a transitive dependency and allow drastically lower performance.

Optionally, the following features can be enabled:

- `bootstrapping` enables gate-bootstrapping, so that the bootstrapping procedure is performed after each gate. A user may wish to have this disabled to use circuit bootstrapping manually.

## Disclaimer

It is a very rough port and not all of the code is idiomatic Rust. Few of the optimizations of the original library are implemented, such as FFT on Lagrange complex polynomials, however may be added in the future.

This is developed as a research project, so I would not recommend to use this for anything more than experimental use.

## License

I am unsure how to license this as it's a port of the TFHE library, however not everything is a direct port and this library has some added features along with some different naming. The original library is licensed under Apache 2.0.
