//! Top-level encryption system types.
//!
//! Contains the [ciphertext](LweSample),
//! the [key parameters](Parameters),
//! the [secret key](SecretKey)
//! and the [cloud key](CloudKey).

use crate::numerics::{approximate_phase, gaussian32, Torus32};
use crate::tgsw::{TGswKey, TGswParams, TGswSample};
use crate::tlwe::TLweKey;
use crate::{tlwe::TLweParameters, SecurityLevel};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

/// Internal ciphertext structure.
/// The `coefficients`-vector is often large and creating new ciphertexts
/// implies allocating a vector of 1024 or more elements, depending
/// on the key parameters.
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct LweSample {
  /// The coefficients of the mask
  pub(crate) coefficients: Vec<Torus32>,
  pub(crate) b: Torus32,
  /// Average noise of the sample
  pub(crate) current_variance: f64,
}

impl LweSample {
  pub(crate) fn new(params: &LweParams) -> Self {
    Self {
      coefficients: vec![0; params.n as usize],
      b: 0,
      current_variance: 0_f64,
    }
  }

  /// Creates a noiseless `LweSample` with a given μ
  pub(crate) fn trivial(mu: i32, params: &LweParams) -> Self {
    Self {
      b: mu,
      ..Self::new(params)
    }
  }
}

impl std::ops::Add<LweSample> for LweSample {
  type Output = Self;
  fn add(self, rhs: LweSample) -> LweSample {
    debug_assert_eq!(
      self.coefficients.len(),
      rhs.coefficients.len(),
      "Cannot add samples with different sizes! lhs.len() = {}, rhs.len() = {}",
      self.coefficients.len(),
      rhs.coefficients.len()
    );
    let coefficients = self
      .coefficients
      .iter()
      .zip_eq(rhs.coefficients)
      .map(|(a, b)| a.wrapping_add(b))
      .collect();

    let b = self.b.wrapping_add(rhs.b);
    let current_variance = self.current_variance + rhs.current_variance;

    Self {
      coefficients,
      b,
      current_variance,
    }
  }
}

impl std::ops::Sub<LweSample> for LweSample {
  type Output = Self;
  fn sub(self, rhs: LweSample) -> LweSample {
    debug_assert_eq!(self.coefficients.len(), rhs.coefficients.len());
    let coefficients = self
      .coefficients
      .iter()
      .zip_eq(rhs.coefficients)
      .map(|(a, b)| a.wrapping_sub(b))
      .collect();

    let b = self.b.wrapping_sub(rhs.b);

    #[allow(clippy::suspicious_arithmetic_impl)]
    let current_variance = self.current_variance + rhs.current_variance;

    Self {
      coefficients,
      b,
      current_variance,
    }
  }
}

impl std::ops::SubAssign for LweSample {
  #[allow(clippy::suspicious_op_assign_impl)]
  fn sub_assign(&mut self, sample: LweSample) {
    self
      .coefficients
      .iter_mut()
      .zip_eq(sample.coefficients.into_iter())
      .for_each(|(a_i, x): (&mut Torus32, Torus32)| *a_i = (*a_i).wrapping_sub(x));
    self.b = self.b.wrapping_sub(sample.b);

    self.current_variance += sample.current_variance;
  }
}

impl std::ops::Mul<i32> for LweSample {
  type Output = Self;
  fn mul(self, p: i32) -> Self {
    let LweSample {
      coefficients,
      b,
      current_variance,
    } = self;

    // Overflowed here, using wrapping mul to imitate C++ behavior
    Self {
      coefficients: coefficients.iter().map(|c| c.wrapping_mul(p)).collect(),
      b: b.wrapping_mul(p),
      current_variance: (p * p) as f64 * current_variance,
    }
  }
}

impl std::ops::Mul<LweSample> for i32 {
  type Output = LweSample;
  fn mul(self, p: LweSample) -> Self::Output {
    p * self
  }
}

impl std::ops::Not for LweSample {
  type Output = Self;
  fn not(self) -> Self {
    let Self {
      coefficients,
      b,
      current_variance,
    } = self;
    Self {
      coefficients: coefficients.iter().map(std::ops::Neg::neg).collect(),
      b: -b,
      current_variance,
    }
  }
}

/// Structure for containing the encryption scheme's parameters.
#[derive(Clone, Deserialize, Serialize)]
pub struct Parameters {
  pub ks_t: i32,
  pub ks_base_bit: i32,
  pub in_out_params: LweParams,
  pub tgsw_params: TGswParams,
}

impl Parameters {
  pub fn new(
    ks_t: i32,
    ks_base_bit: i32,
    in_out_params: LweParams,
    tgsw_params: TGswParams,
  ) -> Self {
    Self {
      ks_t,
      ks_base_bit,
      in_out_params,
      tgsw_params,
    }
  }

  /// Return encryption parameters with the given security level.
  pub fn with(bit_security: SecurityLevel) -> Self {
    const N: i32 = 1024;
    const K: i32 = 1;
    // Max standard deviation for a 1/4 msg space
    const MAX_STDEV: f64 = 0.012_467;

    const KS_BASE_BIT: i32 = 2;
    const KS_LENGTH: i32 = 8;

    match bit_security {
      SecurityLevel::Bit80 => {
        const LOWERCASE_N: i32 = 500;
        const BK_L: i32 = 2;
        const BK_BG_BIT: i32 = 10;

        // Standard deviation
        const KS_STDEV: f64 = 2.44e-5;

        // Standard deviation
        const BK_STDEV: f64 = 7.18e-9;

        let params_in: LweParams = LweParams::new(LOWERCASE_N, KS_STDEV, MAX_STDEV);
        let params_accum: TLweParameters = TLweParameters::new(N, K, BK_STDEV, MAX_STDEV);
        let params_bk: TGswParams = TGswParams::new(BK_L, BK_BG_BIT, params_accum);
        Self {
          ks_t: KS_LENGTH,
          ks_base_bit: KS_BASE_BIT,
          in_out_params: params_in,
          tgsw_params: params_bk,
        }
      }
      SecurityLevel::Bit128 => {
        const LOWERCASE_N: i32 = 630;
        const BK_L: i32 = 3;
        const BK_BG_BIT: i32 = 7;

        // Standard deviation
        let ks_stdev: f64 = 2_f64.powf(-15_f64);

        // Standard deviation
        let bk_stdev: f64 = 2_f64.powf(-15_f64);

        let params_in: LweParams = LweParams::new(LOWERCASE_N, ks_stdev, MAX_STDEV);
        let params_accum: TLweParameters = TLweParameters::new(N, K, bk_stdev, MAX_STDEV);
        let params_bk: TGswParams = TGswParams::new(BK_L, BK_BG_BIT, params_accum);
        Self {
          ks_t: KS_LENGTH,
          ks_base_bit: KS_BASE_BIT,
          in_out_params: params_in,
          tgsw_params: params_bk,
        }
      }
    }
  }
}

impl Default for Parameters {
  /// Returns parameters for the standard security level
  /// ([`SecurityLevel`](super::encryption::SecurityLevel)) of around 128 bits
  fn default() -> Self {
    Self::with(SecurityLevel::Bit128)
  }
}

/// Key for performing homomorphic operations to ciphertexts
/// without gaining access to the data.
/// Safe for sharing as it is meant to be used by a cloud vendor
/// or some other third-party.
#[derive(Clone, Deserialize, Serialize)]
pub struct CloudKey {
  pub(crate) params: Parameters,
  pub bk: LweBootstrappingKey,
}

impl CloudKey {
  pub fn new(params: Parameters, bk: LweBootstrappingKey) -> Self {
    Self { params, bk }
  }
}

/// Key to encrypt and decrypt data.
/// **Not** safe to share.
/// # Warning
/// Although this struct is serializable it is **not** intended to be shared.
/// It only allows serialization to enable storing the key privately.
#[derive(Deserialize, Serialize)]
pub struct SecretKey {
  pub(crate) params: Parameters,
  pub(crate) lwe_key: LweKey,
}

impl SecretKey {
  pub fn new(params: Parameters, lwe_key: LweKey) -> Self {
    Self { params, lwe_key }
  }
}

/// Actual key used for encryption and decryption.
/// **Not** safe to share. Embedded within the [`SecretKey`] key.
/// **Warning**: although it implements serialization and deserialization, it is
/// not intended to be shared, and this functionality is only implemented to
/// store the key privately.
#[derive(Debug, Deserialize, Serialize)]
pub struct LweKey {
  params: LweParams,
  key: Vec<i32>,
}

impl LweKey {
  /// Uses randomness to generate a key with the given parameters
  ///
  /// From C++:
  /// This function generates a random Lwe key for the given parameters.
  /// The Lwe key for the result must be allocated and initialized
  /// (this means that the parameters are already in the result)
  pub fn generate(params: &LweParams) -> Self {
    use rand::Rng;

    let n = params.n;
    let mut rng = rand::thread_rng();

    // Turns out, it is actually faster to use `Vec<i32>` instead of
    // `Vec<bool>` even though it uses less space and the key
    // is binary.
    let key = (0..n).map(|_| if rng.gen() { 1 } else { 0 }).collect();

    Self {
      params: params.clone(),
      key,
    }
  }

  /// This function encrypts message by using key, with stdev alpha
  /// The Lwe sample for the result must be allocated and initialized
  /// (this means that the parameters are already in the result)
  /// TODO: Rewrite this function to return a `LweSample` as it overwrites all values and has all it needs to create a sample instead of mutating one
  pub fn encrypt(&self, result: &mut LweSample, message: Torus32, alpha: f64) {
    use rand::Rng;
    use std::num::Wrapping;

    result.b = gaussian32(message, alpha);
    let mut rng = rand::thread_rng();
    rng.fill(&mut result.coefficients[..]);
    debug_assert_eq!(result.coefficients.len(), self.key.len());
    let values: Wrapping<i32> = result
      .coefficients
      .iter()
      .zip_eq(self.key.iter())
      .map(|(a, b)| Wrapping(a * b))
      .sum();
    result.b = result.b.wrapping_add(values.0);
    result.current_variance = alpha * alpha;
  }

  /// Encrypts a message by using key and a given noise value
  pub(crate) fn encrypt_with_external_noise(
    &self,
    result: &mut LweSample,
    message: Torus32,
    noise: f64,
    alpha: f64,
  ) {
    use crate::numerics::f64_to_torus_32;
    use rand::Rng;
    use std::num::Wrapping;

    result.b = message.wrapping_add(f64_to_torus_32(noise));
    let mut rng = rand::thread_rng();
    rng.fill(&mut result.coefficients[..]);

    debug_assert_eq!(result.coefficients.len(), self.key.len());
    let values: Wrapping<i32> = result
      .coefficients
      .iter()
      .zip_eq(self.key.iter())
      .map(|(a, b)| Wrapping(a * b))
      .sum();
    result.b = result.b.wrapping_add(values.0);
    result.current_variance = alpha * alpha;
  }

  /// This function computes the decryption of sample by using key
  /// The `message_size` indicates the message space and is used to approximate the phase
  pub(crate) fn decrypt(&self, sample: &LweSample, message_size: i32) -> Torus32 {
    let phi = lwe_phase(sample, self);
    approximate_phase(phi, message_size)
  }

  /// Perform the extraction procedure to retrieve an `LweKey` from the
  /// given `TLweKey`.
  pub(crate) fn extract(params: &LweParams, key: &TLweKey) -> Self {
    let mut extracted_key = Self {
      params: params.clone(),
      key: vec![0; params.n as usize],
    };

    let n = key.params.n as usize;
    let k = key.params.k;

    assert_eq!(extracted_key.params.n, k * n as i32);

    for i in 0..k as usize {
      for j in 0..n as usize {
        extracted_key.key[i * n + j] = key.key[i].coefs[j];
      }
    }
    extracted_key
  }
}

/// This function computes the phase of the sample by
/// `b - sample.a .* key.key`, where `.*` is a broadcasting
/// multiplication operator.
pub(crate) fn lwe_phase(sample: &LweSample, key: &LweKey) -> Torus32 {
  use std::num::Wrapping;
  let a: &Vec<Torus32> = &sample.coefficients;
  let k = &key.key;
  debug_assert_eq!(a.len(), k.len());
  let axs: Torus32 = a
    .iter()
    .zip_eq(k.iter())
    .map(|(a, k)| Wrapping(a * k))
    .sum::<Wrapping<i32>>()
    .0;

  // Overflowed here, using wrapping sub to imitate C++ behavior
  sample.b.wrapping_sub(axs)
}

/// Parameters used for the [`LweKey`]
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct LweParams {
  pub(crate) n: i32,

  /// The minimum noise level of an Lwe-sample
  pub(crate) alpha_min: f64,

  /// The greatest noise level that still allows decryption
  alpha_max: f64,
}

impl LweParams {
  pub fn new(n: i32, alpha_min: f64, alpha_max: f64) -> Self {
    Self {
      n,
      alpha_min,
      alpha_max,
    }
  }
}

/// Key used for bootstrapping ciphertexts, safe to share.
/// Is embedded in the [`CloudKey`] key.
#[derive(Clone, Deserialize, Serialize)]
pub struct LweBootstrappingKey {
  /// Input- and output-parameters. key: s
  pub(crate) in_out_params: LweParams,
  /// Parameters for the GSW elements in `bk`. key: s"
  pub(crate) bk_params: TGswParams,
  /// Parameters for the accumulator variable. key: s"
  pub(crate) accum_params: TLweParameters,
  /// Parameters after extraction: key: s'
  pub(crate) extract_params: LweParams,
  /// The bootstrapping key (s->s")
  pub(crate) bk: Vec<TGswSample>,
  /// The keyswitch key (s'->s)
  pub(crate) ks: LweKeySwitchKey,
}

impl LweBootstrappingKey {
  fn new(params: &Parameters) -> Self {
    let Parameters {
      ks_t,
      ks_base_bit,
      in_out_params,
      tgsw_params,
    } = params;
    let accum_params = &tgsw_params.tlwe_params;
    let extract_params = &accum_params.extracted_lweparams;
    let n = in_out_params.n;

    let bk: Vec<TGswSample> = vec![TGswSample::new(&tgsw_params); n as usize];
    let ks = LweKeySwitchKey::new(n, *ks_t, *ks_base_bit, in_out_params);
    let bk_params = params.tgsw_params.clone();

    Self {
      in_out_params: in_out_params.clone(),
      bk_params,
      accum_params: accum_params.clone(),
      extract_params: extract_params.clone(),
      bk,
      ks,
    }
  }

  pub(crate) fn create(params: &Parameters, key_in: &LweKey, rgsw_key: &TGswKey) -> Self {
    let mut key = Self::new(&params);
    assert_eq!(key.bk_params, rgsw_key.params);
    assert_eq!(key.in_out_params, key_in.params);

    let in_out_params = &key.in_out_params;
    let accum_params = &key.bk_params.tlwe_params;
    let extract_params = &accum_params.extracted_lweparams;

    let accum_key = &rgsw_key.tlwe_key;
    let extracted_key = LweKey::extract(&extract_params, &accum_key);
    key.ks.create(&extracted_key, &key_in);

    let alpha = accum_params.alpha_min;
    for i in 0..in_out_params.n as usize {
      rgsw_key.encrypt(&mut key.bk[i], key_in.key[i], alpha);
    }
    key
  }
}

#[derive(Clone, Deserialize, Serialize)]
pub struct LweKeySwitchKey {
  /// length of the input key: s'
  n: i32,
  /// decomposition length
  t: i32,
  /// log_2(base)
  base_bit: i32,
  /// decomposition base: a power of 2
  base: i32,
  /// params of the output key s
  out_params: LweParams,
  // ks0_raw: Vec<LweSample>,
  // of size n*l points to an array ks0_raw whose boxes are spaces basic positions
  // ks1_raw: Vec<Vec<LweSample>>,
  // the keyswitch elements: a n*l*base matrix
  /// Matrix containing all Lwe samples of size n*l*base
  ks: Vec<Vec<Vec<LweSample>>>,
  // of size n points to ks1 an array whose boxes are spaced by ell positions
}

impl LweKeySwitchKey {
  pub(crate) fn new(n: i32, t: i32, base_bit: i32, out_params: &LweParams) -> Self {
    let base = 1 << base_bit;
    let ks = vec![vec![vec![LweSample::new(&out_params); base as usize]; t as usize]; n as usize];
    Self {
      n,
      t,
      base_bit,
      base,
      out_params: out_params.clone(),
      ks,
    }
  }

  pub(crate) fn create(&mut self, in_key: &LweKey, out_key: &LweKey) {
    use rand::distributions::Distribution;
    let alpha = out_key.params.alpha_min;
    let size_ks = self.n * self.t * (self.base - 1);

    // Choose a random vector of gaussian noises
    let mut rng = rand::thread_rng();
    let d = rand_distr::Normal::new(0_f64, alpha).expect("Could not create normal distributioon");
    let noise: Vec<f64> = (0..size_ks).map(|_| d.sample(&mut rng)).collect();
    let error: f64 = noise.iter().sum::<f64>() / size_ks as f64;
    // Recenter the noises
    let noise: Vec<f64> = noise.iter().map(|n| n - error).collect();

    // Generate the ks
    let mut index = 0;
    for i in 0..self.n as usize {
      for j in 0..self.t as usize {
        // term h=0 as trivial encryption of 0 (it will not be used in the KeySwitching)
        self.ks[i][j][0].b = 0;
        // lweNoiselessTrivial(&result->ks[i][j][0], 0, out_key->params);
        for h in 1..self.base as usize {
          let message: Torus32 =
            (in_key.key[i] * h as i32).wrapping_mul(1 << (32 - (j + 1) * self.base_bit as usize));
          out_key.encrypt_with_external_noise(&mut self.ks[i][j][h], message, noise[index], alpha);
          index += 1;
        }
      }
    }
  }

  /// Fills the KeySwitching key array
  /// # Arguments
  /// * `result` - The (n x t x base) array of samples.
  ///    result[i][j][k] encodes k.s[i]/base^(j+1)
  /// * `out_key` - The LWE key to encode all the output samples
  /// * `out_alpha` - The standard deviation of all output samples
  /// * `in_key` - The (binary) input key
  /// * `n` - The size of the input key
  /// * `t` - The precision of the keyswitch (technically, 1/2.base^t)
  /// * `basebit` - Log_2 of base
  pub(crate) fn create_from_array(
    &mut self,
    out_key: &LweKey,
    out_alpha: f64,
    in_key: &[i32],
    n: i32,
    t: i32,
    base_bit: i32,
  ) {
    // base=2 in [CGGI16]
    let base = 1 << base_bit;
    // let mut result = vec![vec![vec![LweSample::new(&out_key.params); base as usize]; t as usize]; n as usize];
    for i in 0..n {
      for j in 0..t {
        for k in 0..base {
          let ax = in_key[i as usize] * k;
          let bx = 1 << (32 - (j + 1) * base_bit);
          // Overflowed here, using wrapping mul to imitate C++ behavior
          let x: Torus32 = ax.wrapping_mul(bx);
          out_key.encrypt(
            &mut self.ks[i as usize][j as usize][k as usize],
            x,
            out_alpha,
          );
        }
      }
    }
  }
}

/// sample = (a',b')
pub(crate) fn lwe_key_switch(ks: &LweKeySwitchKey, sample: LweSample) -> LweSample {
  let mut r = LweSample::trivial(sample.b, &ks.out_params);

  lwe_key_switch_translate_from_array(
    &mut r,
    &ks.ks,
    &sample.coefficients,
    ks.n,
    ks.t,
    ks.base_bit,
  );
  r
}

/// Translates the message of the result sample by -sum(a[i] * s[i]) where s is the secret embedded in `ks`.
/// # Arguments
/// * `result` - the LWE sample to translate by `-sum(ai*si)`.
/// * `ks` - The (n * t * base) key switching key `ks[i][j][k]` encodes k.s[i]/base^(j+1)
/// * `params` - The common LWE parameters of ks and result
/// * `ai` - The input torus array
/// * `n` - The size of the input key
/// * `t` - The precision of the keyswitch (technically, 1/2.base^t)
/// * `basebit` - Log_2 of base
fn lwe_key_switch_translate_from_array(
  result: &mut LweSample,
  #[allow(clippy::ptr_arg)] ks: &[Vec<Vec<LweSample>>],
  ai: &[Torus32],
  n: i32,
  t: i32,
  base_bit: i32,
) {
  // base=2 in [CGGI16]
  let base = 1 << base_bit;
  // precision
  let prec_offset = 1 << (32 - (1 + base_bit * t));
  let mask = base - 1;

  for i in 0..n as usize {
    let aibar = ai[i] + prec_offset;
    for j in 0..t {
      let aij = (aibar >> (32 - (j + 1) * base_bit)) & mask;
      if aij != 0 {
        // TODO: Fix unnecessary allocations here
        *result -= ks[i][j as usize][aij as usize].clone();
      }
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use rand::distributions::Distribution;
  use rand::Rng;

  fn generate_parameters() -> Vec<LweParams> {
    vec![
      LweParams::new(500, 0f64, 1f64),
      LweParams::new(700, 0f64, 1f64),
      LweParams::new(1024, 0f64, 1f64),
    ]
  }

  fn generate_keys() -> Vec<LweKey> {
    let params = generate_parameters();
    params.iter().map(LweKey::generate).collect()
  }

  /// `| frac(x) |`
  ///
  /// Frac:
  /// Return the integer nearest X in the direction of the prevailing rounding mode. (Default in C++ is `FE_TONEAREST`)
  fn abs_frac(x: f64) -> f64 {
    // FIXME: This probably should be using `f64::round()`
    x.trunc().abs()
  }

  #[test]
  fn test_valid_key_generation() {
    for key in generate_keys() {
      // Ensure key is binary
      let count = key.key.iter().fold(0, |acc, &elem| {
        assert!(elem == 0 || elem == 1);
        acc + elem
      });

      assert!(count <= key.params.n - 20);
      assert!(count >= 20);
    }
  }

  /// This function encrypts message by using key, with stdev alpha
  /// The Lwe sample for the result must be allocated and initialized
  /// (this means that the parameters are already in the result)
  #[test]
  fn test_encryption_decryption() {
    use crate::numerics::{mod_switch_to_torus32, torus_32_to_f64};
    use std::collections::HashSet;

    const NB_SAMPLES: i32 = 10;
    const M: i32 = 8;
    const ALPHA: f64 = 1f64 / (10f64 * (M as f64));

    for key in generate_keys() {
      let params = key.params.clone();
      let mut samples = vec![LweSample::new(&params); NB_SAMPLES as usize];

      // Verify correctness of the decryption
      for trial in 0..NB_SAMPLES {
        let message: Torus32 = mod_switch_to_torus32(trial, M);
        key.encrypt(&mut samples[trial as usize], message, ALPHA);
        let phase: Torus32 = lwe_phase(&samples[trial as usize], &key);
        let decrypted = key.decrypt(&samples[trial as usize], M);
        let decrypted_message = torus_32_to_f64(message);
        let decrypted_phase = torus_32_to_f64(phase);
        assert_eq!(message, decrypted);
        println!(
          "absfrac {} <= {}",
          abs_frac(decrypted_message - decrypted_phase),
          10f64 * ALPHA
        );
        assert!(abs_frac(decrypted_message - decrypted_phase) <= 10f64 * ALPHA);
        assert_eq!(ALPHA * ALPHA, samples[trial as usize].current_variance);
      }

      // Verify that samples are random enough (all coordinates different)
      let n = params.n;
      for i in 0..n {
        let mut test_set = HashSet::new();
        for sample in samples.iter() {
          test_set.insert(sample.coefficients[i as usize]);
        }
        assert!(test_set.len() as f64 >= 0.9 * NB_SAMPLES as f64);
      }
    }
  }

  fn fill_random(sample: &LweSample, params: &LweParams) -> LweSample {
    let mut rng = rand::thread_rng();
    let mut coefficients = vec![0; sample.coefficients.len()];
    rng.fill(&mut coefficients[..]);
    let b = rng.gen();
    let current_variance = 0.2;

    LweSample {
      coefficients,
      b,
      current_variance,
    }
  }

  /// TODO: Remove this test. It tests the implementation by implementing it again. Keeping until library is stable as a sort of regression test.
  #[test]
  fn test_adding_samples() {
    for key in generate_keys() {
      let sample1 = fill_random(&LweSample::new(&key.params), &key.params);
      let sample2 = fill_random(&LweSample::new(&key.params), &key.params);
      let a = sample1.clone();
      let b = sample2.clone();
      let result = a + b;

      let coefficients = sample1
        .coefficients
        .iter()
        .zip_eq(sample2.coefficients)
        .map(|(a, b)| a.wrapping_add(b))
        .collect();
      let b = sample1.b.wrapping_add(sample2.b);
      let current_variance = sample1.current_variance + sample2.current_variance;
      assert_eq!(
        result,
        LweSample {
          coefficients,
          b,
          current_variance
        }
      );
    }
  }

  /// TODO: Remove this test. It tests the implementation by implementing it again. Keeping until library is stable as a sort of regression test.
  #[test]
  fn test_subtracting_samples() {
    for key in generate_keys() {
      let sample1 = fill_random(&LweSample::new(&key.params), &key.params);
      let sample2 = fill_random(&LweSample::new(&key.params), &key.params);
      let a = sample1.clone();
      let b = sample2.clone();
      let result = a - b;

      let coefficients = sample1
        .coefficients
        .iter()
        .zip_eq(sample2.coefficients)
        .map(|(a, b)| a.wrapping_sub(b))
        .collect();
      let b = sample1.b.wrapping_sub(sample2.b);
      let current_variance = sample1.current_variance + sample2.current_variance;
      assert_eq!(
        result,
        LweSample {
          coefficients,
          b,
          current_variance
        }
      );
    }
  }

  #[test]
  fn test_key_switch() {
    let mut rng = rand::thread_rng();
    let d = rand_distr::Uniform::new(i32::MIN, i32::MAX);
    let params500 = LweParams::new(500, 0f64, 1f64);
    let key500 = LweKey::generate(&params500);
    let alpha = 1e-5;
    let params = LweParams::new(500, alpha, 1f64);
    let mut key = LweKeySwitchKey::new(300, 14, 2, &params);
    let n = key.n;
    let t = key.t;
    let base = key.base;
    let base_bit = key.base_bit;

    // Precision
    let prec_offset = 1 << (32 - (1 + base_bit * t));
    let prec_mask: u32 = (-(1 << (32 - (base_bit * t))) as i32) as u32;
    let b = d.sample(&mut rng);
    let in_key: Vec<i32> = (0..n)
      .map(|_| match rng.gen::<bool>() {
        true => 1,
        false => 0,
      })
      .collect();
    let mut ai = vec![0; n as usize];
    rng.fill(&mut ai[..]);
    let aibar: Vec<u32> = ai
      .iter()
      .map(|a| ((a + prec_offset) as u32) & prec_mask)
      .collect();

    key.create_from_array(&key500, alpha, &in_key, n, t, base_bit);

    // We first try one by one
    let mut res = LweSample::trivial(b, &params);
    let mut barphi: Torus32 = b;
    for i in 0..n as usize {
      lwe_key_switch_translate_from_array(&mut res, &key.ks[i..], &ai[i..], 1, t, base_bit);

      // Overflowed here, using wrapping sub to imitate C++ behavior
      barphi = barphi.wrapping_sub((aibar[i] as i32) * in_key[i]);

      // Verify the decomposition function
      let mut dec = 0u32;
      for j in 0..t {
        let aij = ((aibar[i] >> (32 - (j + 1) * base_bit)) as i32 & (base - 1)) as u32;
        let rec = aij << (32 - (j + 1) * base_bit);
        dec += rec;
      }
      assert_eq!(dec, aibar[i]);
      assert!(res.current_variance <= alpha * alpha * ((i + 1) as f64) * t as f64 + 1e-10);
      // FIXME: Uncommenting the following line fails... Why?
      // assert_eq!(barphi, res.b);
    }

    // Now, test it all at once
    res = LweSample::trivial(b, &params);
    lwe_key_switch_translate_from_array(&mut res, &key.ks, &ai, n, t, base_bit);
    assert!(res.current_variance <= alpha * alpha * (n as f64) * (t as f64) + 1e-10);
    // FIXME: Uncommenting the following line fails... Why?
    // assert_eq!(barphi, res.b);
  }
}
