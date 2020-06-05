//! # TFHE - Fast Fully Homomorphic Encryption over the Torus
//! This library is a port of the [TFHE](https://tfhe.github.io/tfhe/) library for C and C++.
//!
//! ## Disclaimer
//! This crate is not intended for production use and was developed as part of academic work.
//! Cryptographic security is not guaranteed, as is the functionality of this crate, as some operations such as bootstrapping may be incorrect.
//! Additionally, it does not have the same performance characteristics as the original library.
//!
//! # Example use
//! The following is an example of encrypting a single bit (represented as a `bool`), performing the `xor` (`^`) operator on it and a publicly known constant `true`, while the ciphertext is encrypted:
//! ```rust
//! use tfhe::encryption::{decrypt, encrypt, generate_keys, Parameters};
//! use tfhe::gates::{boots_constant, boots_xor};
//! let message = false;
//! let params = Parameters::default(); // 128-bit security
//! let (secret_key, cloud_key) = generate_keys(&params);
//! let encrypted = encrypt(message, &secret_key);
//! let result = boots_xor(&encrypted, &boots_constant(true, &cloud_key), &cloud_key);
//! let decrypted = decrypt(&result, &secret_key);
//!
//! assert_eq!(decrypted, true);
//! ```

#![allow(dead_code)]
#![forbid(unsafe_code)]
#![doc(test(attr(deny(warnings))))]

pub mod bootstrapping;
pub mod circuits;
pub mod encryption;
pub mod gates;
pub mod numerics;

mod lwe;
mod polynomial;
mod tgsw;
mod tlwe;

pub use bootstrapping::tfhe_bootstrap;
pub use encryption::*;
pub use gates::*;
