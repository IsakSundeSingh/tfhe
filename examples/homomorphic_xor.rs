use tfhe::encryption::{decrypt, encrypt, generate_keys, Parameters};
use tfhe::gates::{constant, xor};

fn main() {
  let message = false;
  let params = Parameters::default();
  let (secret_key, cloud_key) = generate_keys(&params);
  let encrypted = encrypt(message, &secret_key);
  let result = xor(&encrypted, &constant(true, &cloud_key), &cloud_key);
  let decrypted = decrypt(&result, &secret_key);
  assert_eq!(decrypted, true);
}
