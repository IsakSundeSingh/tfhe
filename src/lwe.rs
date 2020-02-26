use crate::tlwe::{TLweParameters, Torus32};
use crate::tsgw::{TGswKey, TGswParams, TGswSample, TGswSampleFFT};

use crate::numerics::gaussian32;

pub struct LweSample {
  /// The coefficients of the mask
  coefficients: Vec<Torus32>,
  b: Torus32,
  /// Average noise of the sample
  current_variance: f64,
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
  params: TFHEGateBootstrappingParameterSet,
  bk: LweBootstrappingKey,
  bk_fft: LweBootstrappingKeyFFT,
}

impl TFHEGateBootstrappingCloudKeySet {
  pub fn new(
    params: TFHEGateBootstrappingParameterSet,
    bk: LweBootstrappingKey,
    bk_fft: LweBootstrappingKeyFFT,
  ) -> Self {
    Self { params, bk, bk_fft }
  }
}

pub struct TFheGateBootstrappingSecretKeySet {
  params: TFHEGateBootstrappingParameterSet,
  lwe_key: LweKey,
  tgsw_key: TGswKey,
  cloud: TFHEGateBootstrappingCloudKeySet,
}
impl TFheGateBootstrappingSecretKeySet {
  pub fn new(
    params: TFHEGateBootstrappingParameterSet,
    bk: LweBootstrappingKey,
    bk_fft: LweBootstrappingKeyFFT,
    lwe_key: LweKey,
    tgsw_key: TGswKey,
  ) -> Self {
    let cloud = TFHEGateBootstrappingCloudKeySet::new(params.clone(), bk, bk_fft);
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
  pub fn generate(params: LweParams) -> Self {
    use rand::distributions::Distribution;

    let n = params.n;
    let d = rand::distributions::Uniform::new_inclusive(0, 1);

    // TODO: Use cryptographically safe RNG
    let mut rng = rand::thread_rng();
    let key = (0..n).map(|_| d.sample(&mut rng)).collect();

    Self { params, key }
  }

  /// This function encrypts message by using key, with stdev alpha
  /// The Lwe sample for the result must be allocated and initialized
  /// (this means that the parameters are already in the result)
  pub fn encrypt(&self, result: &mut LweSample, message: Torus32, alpha: f64) {
    use rand::distributions::Distribution;

    let n = self.params.n;
    result.b = gaussian32(message, alpha);
    let d = rand::distributions::Uniform::new(i32::min_value(), i32::max_value());
    let mut rng = rand::thread_rng();
    for i in 0..n {
      result.coefficients[i as usize] = d.sample(&mut rng);
      result.b += result.coefficients[i as usize] * self.key[i as usize];
    }
    result.current_variance = alpha * alpha;
  }

  /*
   * This function encrypts a message by using key and a given noise value
   */

  pub fn encrypt_with_external_noise(
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
    let d = rand::distributions::Uniform::new(i32::min_value(), i32::max_value());
    let mut rng = rand::thread_rng();
    for i in 0..n {
      result.coefficients[i as usize] = d.sample(&mut rng);
      result.b += result.coefficients[i as usize] * self.key[i as usize];
    }
    result.current_variance = alpha * alpha;
  }
}

//this pub structure contains Lwe parameters
//this pub structure is constant (cannot be modified once initialized):
//the pointer to the param can be passed directly
//to all the Lwe keys that use these params.
#[derive(Clone, Debug, PartialEq)]
pub struct LweParams {
  n: i32,
  /// le plus petit bruit tq sur
  alpha_min: f64,
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
  bk: TGswSample,
  /// the keyswitch key (s'->s)
  ks: LweKeySwitchKey,
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
  /// array which contains all Lwe samples of size nlbase
  ks0_raw: Vec<LweSample>,
  /// of size nl points to an array ks0_raw whose boxes are spaces basic positions
  ks1_raw: Vec<LweSample>,
  /// the keyswitch elements: a n.l.base matrix
  ks: Vec<Vec<LweSample>>,
  // of size n points to ks1 an array whose boxes are spaced by ell positions
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

  #[test]
  fn assert_lwe_key() {
    for param in generate_parameters() {
      let key = LweKey::generate(param);
      let n = key.params.n;
      let s = key.key;
      // Ensure key is binary
      let mut count = 0;
      for i in 0..n {
        assert!(s[i as usize] == 0 || s[i as usize] == 1);
        count += s[i as usize];
      }

      // Sort of useless, isn't it?
      assert!(count <= n - 20);
      assert!(count >= 20);
    }
  }

  #[test]
  #[ignore]
  fn test_encryption_decryption() {
    unimplemented!()
  }
}
