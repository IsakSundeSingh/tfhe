//! # TFHE - Fast Fully Homomorphic Encryption over the Torus
//! This library is a port of the [TFHE](https://tfhe.github.io/tfhe/) library for C and C++.
//!
//! ## Disclaimer
//! This crate is not intended for production use and was developed as part of academic work.
//! Cryptographic security is not guaranteed, as is the functionality of this crate, as some operations such as bootstrapping may be incorrect.
//! Additionally, it does not have the same performance characteristics as the original library.
//!
//! # Example use
//! The following is an example of encrypting a single bit (represented as a `bool`), performing the `xor` (`^`) operator on it and a publicly known constant `true`, while the ciphertext is encrypted.:
//! ```rust
//! use tfhe::bootstrapping::{boots_constant, boots_xor, bootstrapping_parameters, decrypt, encrypt, generate_keys};
//! let message = false;
//! let security = 128;
//! let params = bootstrapping_parameters(security);
//! let (secret_key, cloud_key) = generate_keys(&params);
//! let encrypted = encrypt(message, &secret_key);
//! let result = boots_xor(&encrypted, &boots_constant(true, &cloud_key), &cloud_key);
//! let decrypted = decrypt(&result, &secret_key);
//!
//! assert_eq!(decrypted, true);
//! ```

#![allow(dead_code)]
#![forbid(unsafe_code)]

pub mod bootstrap_internals;
pub mod bootstrapping;
pub mod lwe;
pub mod numerics;
mod polynomial;
pub mod tgsw;
pub mod tlwe;
pub use bootstrap_internals::tfhe_bootstrap;
pub use bootstrapping::*;
