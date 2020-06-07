//! This module contains convenience-circuits to allow the user to more easily
//! build larger circuits manually.
//!
//! It does not guarantee the optimal implementation of such a circuit, so
//! manually construct circuits if you require more optimal circuits.

#[allow(clippy::module_inception)]
pub mod circuits;
pub mod utils;
pub use utils::{as_bits, AsBits};
