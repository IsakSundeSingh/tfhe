use crate::tgsw::{TGswKey, TGswParams, TGswSample, TGswSampleFFT};
use crate::tlwe::TLweKey;
use crate::tlwe::{TLweParameters, Torus32};

use crate::numerics::{approximate_phase, gaussian32};

#[derive(Clone, Debug, PartialEq)]
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
      current_variance: 0f64,
    }
  }
}

impl std::ops::Add<LweSample> for LweSample {
  type Output = Self;
  fn add(self, rhs: LweSample) -> LweSample {
    assert_eq!(
      self.coefficients.len(),
      rhs.coefficients.len(),
      "Cannot add samples with different sizes! lhs.len() = {}, rhs.len() = {}",
      self.coefficients.len(),
      rhs.coefficients.len()
    );
    let coefficients = self
      .coefficients
      .iter()
      .zip(rhs.coefficients)
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
    assert_eq!(self.coefficients.len(), rhs.coefficients.len());
    let coefficients = self
      .coefficients
      .iter()
      .zip(rhs.coefficients)
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

impl std::ops::Mul<i32> for LweSample {
  type Output = Self;
  fn mul(self, p: i32) -> Self {
    let LweSample {
      coefficients,
      b,
      current_variance,
    } = self;

    Self {
      coefficients: coefficients.iter().map(|c| c * p).collect(),
      b: b * p,
      current_variance: (p * p) as f64 * current_variance,
    }
  }
}

impl std::ops::Neg for LweSample {
  type Output = Self;
  fn neg(self) -> Self {
    let Self {
      coefficients,
      b,
      current_variance,
    } = self;
    Self {
      coefficients: coefficients.iter().map(|c| -c).collect(),
      b: -b,
      current_variance,
    }
  }
}

#[derive(Clone)]
pub struct TFHEGateBootstrappingParameterSet {
  pub ks_t: i32,
  pub ks_base_bit: i32,
  pub in_out_params: LweParams,
  pub tgsw_params: TGswParams,
}

impl TFHEGateBootstrappingParameterSet {
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
}

pub struct TFHEGateBootstrappingCloudKeySet {
  pub(crate) params: TFHEGateBootstrappingParameterSet,
  bk: LweBootstrappingKey,
  // bk_fft: LweBootstrappingKeyFFT,
}

impl TFHEGateBootstrappingCloudKeySet {
  pub fn new(
    params: TFHEGateBootstrappingParameterSet,
    bk: LweBootstrappingKey,
    // bk_fft: LweBootstrappingKeyFFT,
  ) -> Self {
    Self { params, bk }
  }
}

pub struct TFheGateBootstrappingSecretKeySet {
  pub(crate) params: TFHEGateBootstrappingParameterSet,
  pub(crate) lwe_key: LweKey,
  tgsw_key: TGswKey,
  pub cloud: TFHEGateBootstrappingCloudKeySet,
}
impl TFheGateBootstrappingSecretKeySet {
  pub fn new(
    params: TFHEGateBootstrappingParameterSet,
    bk: LweBootstrappingKey,
    lwe_key: LweKey,
    tgsw_key: TGswKey,
  ) -> Self {
    let cloud = TFHEGateBootstrappingCloudKeySet::new(params.clone(), bk);
    Self {
      params,
      lwe_key,
      tgsw_key,
      cloud,
    }
  }
}

#[derive(Debug)]
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
    use rand::distributions::Distribution;

    let n = params.n;
    let d = rand::distributions::Uniform::new_inclusive(0, 1);

    // TODO: Use cryptographically safe RNG
    let mut rng = rand::thread_rng();
    let key = (0..n).map(|_| d.sample(&mut rng)).collect();

    Self {
      params: params.clone(),
      key,
    }
  }

  /// This function encrypts message by using key, with stdev alpha
  /// The Lwe sample for the result must be allocated and initialized
  /// (this means that the parameters are already in the result)
  pub fn encrypt(&self, result: &mut LweSample, message: Torus32, alpha: f64) {
    use rand::distributions::Distribution;

    let n = self.params.n;
    result.b = gaussian32(message, alpha);
    let d = rand_distr::Uniform::new(i32::min_value(), i32::max_value());
    let mut rng = rand::thread_rng();
    for i in 0..n {
      result.coefficients[i as usize] = d.sample(&mut rng);

      // Overflowed here, using wrapping add to imitate C++ behavior
      result.b = result
        .b
        .wrapping_add(result.coefficients[i as usize] * self.key[i as usize]);
    }
    result.current_variance = alpha * alpha;
  }

  /*
   * This function encrypts a message by using key and a given noise value
   */

  pub(crate) fn encrypt_with_external_noise(
    &self,
    result: &mut LweSample,
    message: Torus32,
    noise: f64,
    alpha: f64,
  ) {
    use crate::numerics::f64_to_torus_32;
    use rand::distributions::Distribution;

    let n = self.params.n;
    result.b = message + f64_to_torus_32(noise);
    let d = rand_distr::Uniform::new(i32::min_value(), i32::max_value());
    let mut rng = rand::thread_rng();
    for i in 0..n {
      result.coefficients[i as usize] = d.sample(&mut rng);

      // Overflowed here, using wrapping add to imitate C++ behavior
      result.b = result
        .b
        .wrapping_add(result.coefficients[i as usize] * self.key[i as usize]);
    }
    result.current_variance = alpha * alpha;
  }

  /**
   * This function computes the decryption of sample by using key
   * The constant Msize indicates the message space and is used to approximate the phase
   */
  pub(crate) fn decrypt(&self, sample: &LweSample, message_size: i32) -> Torus32 {
    let phi = lwe_phase(sample, self);
    approximate_phase(phi, message_size)
  }

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

/**
 * This function computes the phase of sample by using key : phi = b - a.s
 */
pub(crate) fn lwe_phase(sample: &LweSample, key: &LweKey) -> Torus32 {
  let n = key.params.n;
  let mut axs: Torus32 = 0;
  let a: &Vec<Torus32> = &sample.coefficients;
  let k = &key.key;

  for i in 0..n {
    // Overflowed here, using wrapping add to imitate C++ behavior
    axs = axs.wrapping_add(a[i as usize] * k[i as usize]);
  }

  // Overflowed here, using wrapping sub to imitate C++ behavior
  sample.b.wrapping_sub(axs)
}

//this pub structure contains Lwe parameters
//this pub structure is constant (cannot be modified once initialized):
//the pointer to the param can be passed directly
//to all the Lwe keys that use these params.
#[derive(Clone, Debug, PartialEq)]
pub struct LweParams {
  pub(crate) n: i32,
  /// le plus petit bruit tq sur
  pub(crate) alpha_min: f64,
  /// le plus gd bruit qui permet le déchiffrement
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

pub struct LweBootstrappingKey {
  /// paramètre de l'input et de l'output. key: s
  in_out_params: LweParams,
  /// params of the Gsw elems in bk. key: s"
  bk_params: TGswParams,
  /// params of the accum variable key: s"
  accum_params: TLweParameters,
  /// params after extraction: key: s'
  extract_params: LweParams,
  /// the bootstrapping key (s->s")
  bk: Vec<TGswSample>,
  /// the keyswitch key (s'->s)
  ks: LweKeySwitchKey,
}

impl LweBootstrappingKey {
  fn new(params: &TFHEGateBootstrappingParameterSet) -> Self {
    let TFHEGateBootstrappingParameterSet {
      ks_t,
      ks_base_bit,
      in_out_params,
      tgsw_params,
    } = params;
    let accum_params = &tgsw_params.tlwe_params;
    let extract_params = &accum_params.extracted_lweparams;
    let n = in_out_params.n;
    let big_n = extract_params.n;

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

  pub(crate) fn create(
    params: &TFHEGateBootstrappingParameterSet,
    key_in: &LweKey,
    rgsw_key: &TGswKey,
  ) -> Self {
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

pub struct LweBootstrappingKeyFFT {
  ///< paramètre de l'input et de l'output. key: s
  in_out_params: LweParams,
  ///< params of the Gsw elems in bk. key: s"
  bk_params: TGswParams,
  ///< params of the accum variable key: s"
  accum_params: TLweParameters,
  ///< params after extraction: key: s'
  extract_params: LweParams,
  ///< the bootstrapping key (s->s")
  bk_fft: TGswSampleFFT,
  ///< the keyswitch key (s'->s)
  ks: LweKeySwitchKey,
}

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
  /// array which contains all Lwe samples of size n*l*base
  // ks0_raw: Vec<LweSample>,
  // of size n*l points to an array ks0_raw whose boxes are spaces basic positions
  // ks1_raw: Vec<Vec<LweSample>>,
  // the keyswitch elements: a n.l.base matrix
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
    let d = rand_distr::Normal::new(0f64, alpha).expect("Could not create normal distributioon");
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
          //    Torus32 mess = (in_key->key[i] * h)*(1<<(32-(j+1)*basebit));
          let message: Torus32 =
            (in_key.key[i] * h as i32) * (1 << (32 - (j + 1) * self.base_bit as usize));
          out_key.encrypt_with_external_noise(&mut self.ks[i][j][h], message, noise[index], alpha);
          //    lweSymEncryptWithExternalNoise(&result->ks[i][j][h], mess, noise[index], alpha, out_key);
          index += 1;
        }
      }
    }
  }
}

#[cfg(test)]
mod tests {
  use super::*;

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
    // TODO: Figure out if this really should be truncating or using `f64::round()`
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
    use rand::distributions::Distribution;
    let d = rand_distr::Uniform::new(i32::min_value(), i32::max_value());
    let mut rng = rand::thread_rng();
    let coefficients = (0..sample.coefficients.len())
      .map(|_| d.sample(&mut rng))
      .collect();
    let b = d.sample(&mut rng);
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
        .zip(sample2.coefficients)
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
        .zip(sample2.coefficients)
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
}
