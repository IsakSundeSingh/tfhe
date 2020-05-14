use tfhe::encryption::{decrypt, encrypt, generate_keys, Parameters};
use tfhe::gates::{boots_constant, boots_xor};

fn main() {
  let message = false;
  let params = Parameters::default();
  let (secret_key, cloud_key) = generate_keys(&params);
  let encrypted = encrypt(message, &secret_key);
  let result = boots_xor(&encrypted, &boots_constant(true, &cloud_key), &cloud_key);
  let decrypted = decrypt(&result, &secret_key);
  assert_eq!(decrypted, true);
}
