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

## Testing

Tests can be run using `cargo test`.

## Benchmarking

Benchmarks can be run using `cargo bench`.

Some benchmarks require specific features to be enabled. Required features are found in [`Cargo.toml`](Cargo.toml) under each specific benchmark's section, `[[bench]]`.

## Disclaimer

It is a very rough port and not all of the code is idiomatic Rust. Few of the optimizations of the original library are implemented, such as FFT on Lagrange complex polynomials, however may be added in the future.

This is developed as a research project, so I would not recommend to use this for anything more than experimental use.

I am not a cryptologist.

## Contributions

Contributions both in terms of issues and pull requests are welcome.
Please use the issue templates when filling out an issue.
Do note that the project is experimental and has the scars from this as well.

## License

This project is a Rust port of the C++ library [TFHE](https://github.com/tfhe/tfhe), and is licensed under the Apache 2.0 License.
