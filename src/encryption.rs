//! Contains functions and types to generate keys, encrypt and decrypt data.
pub use crate::lwe::{CloudKey, LweSample, Parameters, SecretKey};

use crate::lwe::{LweBootstrappingKey, LweKey};
use crate::{
  numerics::{encode_message, Torus32},
  tgsw::TGswKey,
};

/// Generate default gate bootstrapping parameters.
/// `minimum_lambda` determines the minimum security level
pub fn bootstrapping_parameters(minimum_lambda: i32) -> Parameters {
  if minimum_lambda > 128 {
    panic!("Sorry, for now, the parameters are only implemented for about 128bit of security!");
  }
  Parameters::default()
}

/// Generate a keypair
/// Returns a tuple of the secret key and cloud key.
/// # Example
///
/// ```rust
/// # use tfhe::encryption::{generate_keys, Parameters};
/// let params = Parameters::default();
/// # #[allow(unused)]
/// let (secret_key, cloud_key) = generate_keys(&params);
/// ```
pub fn generate_keys(params: &Parameters) -> (SecretKey, CloudKey) {
  let lwe_key: LweKey = LweKey::generate(&params.in_out_params);
  let tgsw_key = TGswKey::generate(&params.tgsw_params);
  let bk = LweBootstrappingKey::create(&params, &lwe_key, &tgsw_key);
  (
    SecretKey::new(params.clone(), lwe_key),
    CloudKey::new(params.clone(), bk),
  )
}

/// Encrypts a given message with the secret key.
pub fn encrypt(message: bool, key: &SecretKey) -> LweSample {
  const _1S8: Torus32 = encode_message(1, 8);
  let mu: Torus32 = if message { _1S8 } else { -_1S8 };
  let alpha = key.params.in_out_params.alpha_min;
  let mut sample = LweSample::new(&key.params.in_out_params);
  key.lwe_key.encrypt(&mut sample, mu, alpha);
  sample
}

/// Decrypts a given sample into a boolean value using the given secret key.
pub fn decrypt(sample: &LweSample, key: &SecretKey) -> bool {
  crate::lwe::lwe_phase(sample, &key.lwe_key) > 0
}
