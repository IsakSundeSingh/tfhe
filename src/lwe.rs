use crate::tlwe::{TLweParameters, Torus32};
use crate::tsgw::{TGswKey, TGswParams, TGswSample, TGswSampleFFT};

pub struct LweSample {
  /// The coefficients of the mask
  coefficients: Vec<Torus32>,
  b: Torus32,
  /// Average noise of the sample
  current_variance: f64,
}

pub struct TFHEGateBootstrappingParameterSet {
  ks_t: i32,
  ks_base_bit: i32,
  in_out_params: LweParams,
  tgsw_params: TGswParams,
}

pub struct TFHEGateBootstrappingCloudKeySet {
  params: TFHEGateBootstrappingParameterSet,
  bk: LweBootstrappingKey,
  bk_fft: LweBootstrappingKeyFFT,
}

pub struct TFheGateBootstrappingSecretKeySet {
  params: TFHEGateBootstrappingParameterSet,
  lwe_key: LweKey,
  tgsw_key: TGswKey,
  cloud: TFHEGateBootstrappingCloudKeySet,
}

pub struct LweKey {
  params: LweParams,
  key: i32,
  // const LweParams* params;
  // int32_t* key;
}
//this pub structure contains Lwe parameters
//this pub structure is constant (cannot be modified once initialized):
//the pointer to the param can be passed directly
//to all the Lwe keys that use these params.
pub struct LweParams {
  n: i32,
  /// le plus petit bruit tq sur
  alpha_min: f64,
  /// le plus gd bruit qui permet le déchiffrement
  alpha_max: f64,
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
