use tfhe::encryption::{
  decrypt, encrypt, generate_keys, generate_parameters, Parameters, SecurityLevel,
};
use tfhe::gates::{
  boots_and, boots_andny, boots_andyn, boots_constant, boots_nand, boots_nor, boots_not, boots_or,
  boots_orny, boots_oryn, boots_xnor, boots_xor,
};

#[test]
fn test_encrypt_decrypt() {
  let params = Parameters::default();
  let (secret_key, _cloud_key) = generate_keys(&params);
  let every_combo = vec![true, false];

  for x in every_combo {
    let encrypted = encrypt(x, &secret_key);
    let decrypted = decrypt(&encrypted, &secret_key);

    assert_eq!(x, decrypted);
  }
}

#[test]
fn test_bootstrapping_constant() {
  let params = generate_parameters(SecurityLevel::Bit80);
  let (secret_key, cloud_key) = generate_keys(&params);
  let every_combo = vec![true, false];

  for x in every_combo {
    let bootstrapped = boots_constant(x, &cloud_key);
    let decrypted = decrypt(&bootstrapped, &secret_key);
    assert_eq!(x, decrypted);
  }
}

macro_rules! test_binary_gate {
  ($test_name: ident, $binary_gate:ident, $encrypted_gate:ident) => {
    #[test]
    fn $test_name() {
      let params = generate_parameters(SecurityLevel::Bit80);
      let (secret_key, cloud_key) = generate_keys(&params);
      let enc_true = encrypt(true, &secret_key);
      let enc_false = encrypt(false, &secret_key);

      let every_combo = vec![(true, true), (true, false), (false, true), (false, false)];

      for (a, b) in every_combo {
        let enc_a = match a {
          true => &enc_true,
          false => &enc_false,
        };
        let enc_b = match b {
          true => &enc_true,
          false => &enc_false,
        };

        let encrypted = $encrypted_gate(enc_a, enc_b, &cloud_key);
        let decrypted = decrypt(&encrypted, &secret_key);
        assert_eq!(decrypted, $binary_gate(a, b));
      }
    }
  };
}

macro_rules! test_unary_gate {
  ($test_name: ident, $unary_gate:ident, $encrypted_gate:ident) => {
    #[test]
    fn $test_name() {
      let params = generate_parameters(SecurityLevel::Bit80);
      let (secret_key, cloud_key) = generate_keys(&params);
      let enc_true = encrypt(true, &secret_key);
      let enc_false = encrypt(false, &secret_key);

      let every_combo = vec![true, false];

      for x in every_combo {
        let encrypted = match x {
          true => &enc_true,
          false => &enc_false,
        };

        let result = $encrypted_gate(encrypted, &cloud_key);
        let decrypted = decrypt(&result, &secret_key);

        assert_eq!(decrypted, $unary_gate(x));
      }
    }
  };
}

fn and(a: bool, b: bool) -> bool {
  a && b
}

fn or(a: bool, b: bool) -> bool {
  a || b
}

fn nand(a: bool, b: bool) -> bool {
  !(a && b)
}

fn not(x: bool) -> bool {
  !x
}

fn nor(a: bool, b: bool) -> bool {
  !(a || b)
}

fn xor(a: bool, b: bool) -> bool {
  a ^ b
}

fn xnor(a: bool, b: bool) -> bool {
  a == b
}

fn andny(a: bool, b: bool) -> bool {
  (!a) && b
}

fn andyn(a: bool, b: bool) -> bool {
  a && (!b)
}

fn orny(a: bool, b: bool) -> bool {
  (!a) || b
}

fn oryn(a: bool, b: bool) -> bool {
  a || (!b)
}

test_binary_gate!(test_and_gate, and, boots_and);
test_binary_gate!(test_or_gate, or, boots_or);
test_binary_gate!(test_nand_gate, nand, boots_nand);
test_binary_gate!(test_nor_gate, nor, boots_nor);
test_binary_gate!(test_xor_gate, xor, boots_xor);
test_binary_gate!(test_xnor_gate, xnor, boots_xnor);
test_binary_gate!(test_andny_gate, andny, boots_andny);
test_binary_gate!(test_andyn_gate, andyn, boots_andyn);
test_binary_gate!(test_orny_gate, orny, boots_orny);
test_binary_gate!(test_oryn_gate, oryn, boots_oryn);
test_unary_gate!(test_not_gate, not, boots_not);
