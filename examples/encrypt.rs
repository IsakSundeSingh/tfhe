use tfhe::encryption::{decrypt, encrypt, generate_keys, Parameters};

fn main() {
  let message = true;
  let params = Parameters::default();
  let (secret_key, _cloud_key) = generate_keys(&params);
  let encrypted = encrypt(message, &secret_key);
  let decrypted = decrypt(&encrypted, &secret_key);

  assert_eq!(message, decrypted);
}
