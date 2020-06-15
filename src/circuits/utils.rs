//! A small module containing useful utility functions.

/// A helper function for converting a byte to an array of bits.
pub const fn as_bits(byte: u8) -> [bool; 8] {
  let bit7: bool = ((byte >> 7) & 0b1) == 1;
  let bit6: bool = ((byte >> 6) & 0b1) == 1;
  let bit5: bool = ((byte >> 5) & 0b1) == 1;
  let bit4: bool = ((byte >> 4) & 0b1) == 1;
  let bit3: bool = ((byte >> 3) & 0b1) == 1;
  let bit2: bool = ((byte >> 2) & 0b1) == 1;
  let bit1: bool = ((byte >> 1) & 0b1) == 1;
  let bit0: bool = (byte & 0b1) == 1;

  [bit7, bit6, bit5, bit4, bit3, bit2, bit1, bit0]
}

/// Trait that allows representing a byte as an array of bits.
///
/// Const generics blocks this from being a more useful trait that can return
/// `N` bits instead of only `8`.
/// See issue #44580 https://github.com/rust-lang/rust/issues/44580
pub trait AsBits {
  /// Represent the bits of the byte as an array of boolean values (bits).
  /// Array is in big-endian order, where MSB is the first value of the array.
  fn as_bits(self) -> [bool; 8];
}

impl AsBits for u8 {
  fn as_bits(self) -> [bool; 8] {
    as_bits(self)
  }
}
