use tfhe::bootstrapping::{
  boots_sym_decrypt, boots_sym_encrypt, new_default_gate_bootstrapping_parameters,
  new_random_gate_bootstrapping_secret_keyset,
};

fn main() {
  let message = true;
  let security = 128;
  let params = new_default_gate_bootstrapping_parameters(security);
  let secret_key = new_random_gate_bootstrapping_secret_keyset(&params);
  let encrypted = boots_sym_encrypt(message, &secret_key);
  let decrypted = boots_sym_decrypt(&encrypted, &secret_key);

  assert_eq!(message, decrypted);
}
