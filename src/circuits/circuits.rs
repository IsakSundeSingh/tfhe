//! Contains commonly-used circuits.
//!
//! # Note
//! This module is solely for convenience usage. It does not implement
//! efficient circuits of any kind, and many of the circuits listed are
//! possible to implement in much more optimal ways.
//! The user is therefore encouraged to implement the required circuits
//! themselves, as they have better control of the use-cases and can
//! tailor the circuits to their needs.

use crate::{and, constant, mux, or, xor, CloudKey, LweSample};
use itertools::Itertools;

/// Comparison gate.
/// Compares two ciphertexts using greater than or equal to-comparison.
/// Effectively does `a <= b`
/// # Panics
/// Panics if `a` and `b` are not the same length
/// # Example
/// ```no_run
/// # use tfhe::circuits::compare;
/// # use tfhe::encryption::{encrypt, decrypt, generate_keys, Parameters};
/// # let params = Parameters::default();
/// # let (secret_key, cloud_key) = generate_keys(&params);
/// let a = vec![false, false];
/// let b = vec![false, true]; // b is greater than a
/// let encrypted_a: Vec<_> = a.into_iter().map(|x| encrypt(x, &secret_key)).collect();
/// let encrypted_b: Vec<_> = b.into_iter().map(|x| encrypt(x, &secret_key)).collect();
/// let encrypted_comparison = compare(&encrypted_a[..], &encrypted_b[..], &cloud_key);
/// let comparison = decrypt(&encrypted_comparison, &secret_key);
/// assert_eq!(comparison, true); // `a <= b` is true
/// ```
pub fn compare(a: &[LweSample], b: &[LweSample], key: &CloudKey) -> LweSample {
  let one = constant(true, key);
  let mut carry = one;
  assert_eq!(a.len(), b.len());
  for i in 0..a.len() {
    carry = le(&a[i], &b[i], &carry, key);
  }
  carry
}

/// Less than or equal to comparator with carry input.
fn le(a: &LweSample, b: &LweSample, rim: &LweSample, key: &CloudKey) -> LweSample {
  let condition = xor(a, b, key);
  mux(&condition, b, rim, key)
}

/// Homomorphic bitwise-equality operator
pub fn eq(a: &[LweSample], b: &[LweSample], key: &CloudKey) -> LweSample {
  assert!(!a.is_empty());
  assert_eq!(a.len(), b.len());

  let mut bits = vec![];
  for (bit_a, bit_b) in a.iter().zip_eq(b.iter()) {
    bits.push(eq_bit(bit_a, bit_b, key));
  }

  bits
    .iter()
    .fold(constant(true, key), |res, bit| and(&res, bit, key))
}

/// Homomorphic bit-equality operator
pub fn eq_bit(a: &LweSample, b: &LweSample, key: &CloudKey) -> LweSample {
  let res = xor(a, b, key);
  let one = constant(true, key);
  xor(&res, &one, key)
}

pub fn swap(_a: &LweSample, _b: &LweSample, _key: &CloudKey) -> LweSample {
  todo!()
}

/// Binary adder accepting two single encrypted binary digits.
/// Returns the sum and the carry, respectively.
pub fn half_adder(a: &LweSample, b: &LweSample, key: &CloudKey) -> (LweSample, LweSample) {
  (xor(a, b, key), and(a, b, key))
}

/// Binary full adder accepting two input signals and the initial carry-input.
/// Returns the sum and the carry
pub fn full_adder(
  a: &LweSample,
  b: &LweSample,
  carry_in: &LweSample,
  key: &CloudKey,
) -> (LweSample, LweSample) {
  let a_b_xor = xor(a, b, key);
  let sum = xor(&a_b_xor, carry_in, key);
  let carry = {
    let carry_in_and_a_xor_b = and(carry_in, &a_b_xor, key);
    let a_and_b = and(a, b, key);
    or(&carry_in_and_a_xor_b, &a_and_b, key)
  };
  (sum, carry)
}

/// Ripple-carry adder that accepts two n-bit encrypted numbers and computes
/// the encrypted sum.
/// Returns the n-bit sum along with the carry bit.
///
/// Example diagram explaining this circuit is found at:
/// https://en.wikipedia.org/wiki/Adder_(electronics)#/media/File:4-bit_ripple_carry_adder.svg
pub fn add(a: &[LweSample], b: &[LweSample], key: &CloudKey) -> (Vec<LweSample>, LweSample) {
  assert_eq!(
    a.len(),
    b.len(),
    "Cannot add two numbers with different number of bits!"
  );

  a.iter().zip_eq(b.iter()).fold(
    (Vec::with_capacity(a.len()), constant(false, key)),
    |(mut sum_bits, carry_in), (a, b)| {
      let (sum, carry_out) = full_adder(a, b, &carry_in, key);
      sum_bits.push(sum);
      (sum_bits, carry_out)
    },
  )
}
