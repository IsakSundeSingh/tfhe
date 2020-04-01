use tfhe::bootstrapping::{bootstrapping_parameters, decrypt, encrypt, generate_keys};

fn main() {
  let message = true;
  let security = 128;
  let params = bootstrapping_parameters(security);
  let (secret_key, _cloud_key) = generate_keys(&params);
  let encrypted = encrypt(message, &secret_key);
  let decrypted = decrypt(&encrypted, &secret_key);

  assert_eq!(message, decrypted);
}
