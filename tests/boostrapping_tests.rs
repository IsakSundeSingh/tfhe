use tfhe::bootstrapping::{
  boots_sym_decrypt, boots_sym_encrypt, new_default_gate_bootstrapping_parameters,
  new_gate_bootstrapping_ciphertext, new_random_gate_bootstrapping_secret_keyset,
};

#[test]
fn test_encrypt_decrypt_true() {
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
