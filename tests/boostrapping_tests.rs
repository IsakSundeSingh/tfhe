use tfhe::bootstrapping::boots_xnor;
use tfhe::bootstrapping::{
  boots_and, boots_nand, boots_nor, boots_not, boots_or, boots_sym_decrypt, boots_sym_encrypt,
  boots_xor, new_default_gate_bootstrapping_parameters,
  new_random_gate_bootstrapping_secret_keyset,
};

#[test]
fn test_encrypt_decrypt_true_is_true() {
  let message = true;
  let security = 128;
  let params = new_default_gate_bootstrapping_parameters(security);
  let secret_key = new_random_gate_bootstrapping_secret_keyset(&params);
  let encrypted = boots_sym_encrypt(message, &secret_key);
  let decrypted = boots_sym_decrypt(&encrypted, &secret_key);

  assert_eq!(message, decrypted);
}

#[test]
fn test_encrypt_decrypt_false_is_false() {
  let message = false;
  let security = 128;
  let params = new_default_gate_bootstrapping_parameters(security);
  let secret_key = new_random_gate_bootstrapping_secret_keyset(&params);
  let encrypted = boots_sym_encrypt(message, &secret_key);
  let decrypted = boots_sym_decrypt(&encrypted, &secret_key);

  assert_eq!(message, decrypted);
}

macro_rules! test_binary_gate {
  ($test_name: ident, $binary_gate:ident, $encrypted_gate:ident) => {
    #[test]
    fn $test_name() {
      let security = 128;
      let params = new_default_gate_bootstrapping_parameters(security);
      let secret_key = new_random_gate_bootstrapping_secret_keyset(&params);
      let cloud_key = &secret_key.cloud;
      let enc_true = boots_sym_encrypt(true, &secret_key);
      let enc_false = boots_sym_encrypt(false, &secret_key);

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
        let decrypted = boots_sym_decrypt(&encrypted, &secret_key);
        assert_eq!(decrypted, $binary_gate(a, b));
      }
    }
  };
}

macro_rules! test_unary_gate {
  ($test_name: ident, $unary_gate:ident, $encrypted_gate:ident) => {
    #[test]
    fn $test_name() {
      let security = 128;
      let params = new_default_gate_bootstrapping_parameters(security);
      let secret_key = new_random_gate_bootstrapping_secret_keyset(&params);
      let cloud_key = &secret_key.cloud;
      let enc_true = boots_sym_encrypt(true, &secret_key);
      let enc_false = boots_sym_encrypt(false, &secret_key);

      let every_combo = vec![true, false];

      for x in every_combo {
        let encrypted = match x {
          true => &enc_true,
          false => &enc_false,
        };

        let result = $encrypted_gate(encrypted, &cloud_key);
        let decrypted = boots_sym_decrypt(&result, &secret_key);

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

test_binary_gate!(test_and_gate, and, boots_and);
test_binary_gate!(test_or_gate, or, boots_or);
test_binary_gate!(test_nand_gate, nand, boots_nand);
test_binary_gate!(test_nor_gate, nor, boots_nor);
test_binary_gate!(test_xor_gate, xor, boots_xor);
test_binary_gate!(test_xnor_gate, xnor, boots_xnor);
test_unary_gate!(test_not_gate, not, boots_not);
