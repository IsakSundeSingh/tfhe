//! Contains commonly-used circuits.

use crate::{boots_and, boots_constant, boots_mux, boots_xor, CloudKey, LweSample};

/// Comparison gate.
/// Compares two ciphertexts using greater than or equal to-comparison.
/// Effectively does `a <= b`
/// # Panics
/// Panics if `a` and `b` are not the same length
pub fn compare(a: &[LweSample], b: &[LweSample], key: &CloudKey) -> LweSample {
  let one = boots_constant(true, key);
  let mut carry = one;
  assert_eq!(a.len(), b.len());
  for i in 0..a.len() {
    carry = le(&a[i], &b[i], &carry, key);
  }
  carry
}

/// Less than or equal to comparator with carry input.
fn le(a: &LweSample, b: &LweSample, rim: &LweSample, key: &CloudKey) -> LweSample {
  let condition = boots_xor(a, b, key);
  boots_mux(&condition, b, rim, key)
}

/// Homomorphic bitwise-equality operator
pub fn eq(a: &[LweSample], b: &[LweSample], key: &CloudKey) -> LweSample {
  assert!(!a.is_empty());
  assert_eq!(a.len(), b.len());

  let mut bits = vec![];
  for (bit_a, bit_b) in a.iter().zip(b.iter()) {
    bits.push(eq_bit(bit_a, bit_b, key));
  }

  bits.iter().fold(boots_constant(true, key), |res, bit| {
    boots_and(&res, bit, key)
  })
}

/// Homomorphic bit-equality operator
pub fn eq_bit(a: &LweSample, b: &LweSample, key: &CloudKey) -> LweSample {
  let res = boots_xor(a, b, key);
  let one = boots_constant(true, key);
  boots_xor(&res, &one, key)
}

pub fn swap(_a: &LweSample, _b: &LweSample, _key: &CloudKey) -> LweSample {
  todo!()
}
