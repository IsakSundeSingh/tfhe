use crate::lwe::{
  LweBootstrappingKey, LweKey, LweParams, LweSample, TFHEGateBootstrappingCloudKeySet,
  TFHEGateBootstrappingParameterSet, TFheGateBootstrappingSecretKeySet,
};
use crate::tgsw::{TGswKey, TGswParams};
use crate::tlwe::TLweParameters;

use crate::numerics::mod_switch_to_torus32;
//////////////////////////////////////////
// Gate bootstrapping public interface
//////////////////////////////////////////
/** generate default gate bootstrapping parameters */
pub fn new_default_gate_bootstrapping_parameters(
  minimum_lambda: i32,
) -> TFHEGateBootstrappingParameterSet {
  if minimum_lambda > 128 {
    panic!("Sorry, for now, the parameters are only implemented for about 128bit of security!");
  }

  const N: i32 = 1024;
  const K: i32 = 1;
  const LOWERCASE_N: i32 = 500;
  const BK_L: i32 = 2;
  const BK_BG_BIT: i32 = 10;
  const KS_BASE_BIT: i32 = 2;
  const KS_LENGTH: i32 = 8;

  // Standard deviation
  const KS_STDEV: f64 = 2.44e-5;

  // Standard deviation
  const BK_STDEV: f64 = 7.18e-9;

  // Max standard deviation for a 1/4 msg space
  const MAX_STDEV: f64 = 0.012_467;

  let params_in: LweParams = LweParams::new(LOWERCASE_N, KS_STDEV, MAX_STDEV);
  let params_accum: TLweParameters = TLweParameters::new(N, K, BK_STDEV, MAX_STDEV);
  let params_bk: TGswParams = TGswParams::new(BK_L, BK_BG_BIT, params_accum);

  TFHEGateBootstrappingParameterSet::new(KS_LENGTH, KS_BASE_BIT, params_in, params_bk)
}

/** generate a random gate bootstrapping secret key */
pub fn new_random_gate_bootstrapping_secret_keyset(
  params: &TFHEGateBootstrappingParameterSet,
) -> TFheGateBootstrappingSecretKeySet {
  let lwe_key: LweKey = LweKey::generate(&params.in_out_params);
  let tgsw_key = TGswKey::generate(&params.tgsw_params);
  let bk = LweBootstrappingKey::create(&params, &lwe_key, &tgsw_key);
  TFheGateBootstrappingSecretKeySet::new(params.clone(), bk, lwe_key, tgsw_key)
}

/** generate a new unititialized ciphertext */
pub fn new_gate_bootstrapping_ciphertext(params: &TFHEGateBootstrappingParameterSet) -> LweSample {
  LweSample::new(&params.in_out_params)
}

/** generate a new unititialized ciphertext array of length nbelems */
pub fn new_gate_bootstrapping_ciphertext_array(
  nbelems: i32,
  params: &TFHEGateBootstrappingParameterSet,
) -> Vec<LweSample> {
  vec![new_gate_bootstrapping_ciphertext(&params); nbelems as usize]
}

/** encrypts a boolean */
pub fn boots_sym_encrypt(message: bool, key: &TFheGateBootstrappingSecretKeySet) -> LweSample {
  let _1s8 = crate::numerics::mod_switch_to_torus32(1, 8);
  let mu = if message { _1s8 } else { -_1s8 };
  let alpha = key.params.in_out_params.alpha_min;
  let mut sample = LweSample::new(&key.params.in_out_params);
  key.lwe_key.encrypt(&mut sample, mu, alpha);
  sample
}

/** decrypts a boolean */
pub fn boots_sym_decrypt(sample: &LweSample, key: &TFheGateBootstrappingSecretKeySet) -> bool {
  use crate::lwe::lwe_phase;
  lwe_phase(sample, &key.lwe_key) > 0
}

/** bootstrapped Constant (true or false) trivial Gate */
pub fn boots_constant(value: i32, bk: &TFHEGateBootstrappingCloudKeySet) -> LweSample {
  unimplemented!()
}

/** bootstrapped Nand Gate */
pub fn boots_nand(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  let mu = mod_switch_to_torus32(1, 8);
  let in_out_params = &bk.params.in_out_params;

  // Compute: (0,1/8) - ca - cb
  let nand = mod_switch_to_torus32(1, 8);
  let temp_result = LweSample {
    coefficients: vec![0; in_out_params.n as usize],
    b: nand,
    current_variance: 0f64,
  };

  temp_result - ca.clone() - cb.clone()

  // If the phase is positive, the result is 1/8,
  // otherwise the result is -1/8
  // TODO: Actually implement the bootstrapping so gates can be chained!
  // tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
}

/** bootstrapped Or Gate:  */
pub fn boots_or(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  let mu = mod_switch_to_torus32(1, 8);
  let in_out_params = &bk.params.in_out_params;

  // Compute: (0,1/8) + ca + cb
  let or = mod_switch_to_torus32(1, 8);

  let temp_result = LweSample {
    coefficients: vec![0; in_out_params.n as usize],
    b: or,
    current_variance: 0f64,
  };

  temp_result + ca.clone() + cb.clone()

  // If the phase is positive, the result is 1/8,
  // otherwise the result is -1/8
  // TODO: Actually implement the bootstrapping so gates can be chained!
  // tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
}

/** bootstrapped And Gate: result = a and b */
pub fn boots_and(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  let mu = mod_switch_to_torus32(1, 8);
  let in_out_params = &bk.params.in_out_params;

  // Compute: (0,-1/8) + ca + cb
  let and = mod_switch_to_torus32(-1, 8);

  let temp_result = LweSample {
    coefficients: vec![0; in_out_params.n as usize],
    b: and,
    current_variance: 0f64,
  };

  temp_result + ca.clone() + cb.clone()

  // If the phase is positive, the result is 1/8,
  // otherwise the result is -1/8
  // TODO: Actually implement the bootstrapping so gates can be chained!
  // tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
}

/** bootstrapped Xor Gate: result = a xor b */
pub fn boots_xor(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  let mu = mod_switch_to_torus32(1, 8);
  let in_out_params = &bk.params.in_out_params;

  // Compute: (0,1/4) + 2*(ca + cb)
  let xor = mod_switch_to_torus32(1, 4);

  let temp_result = LweSample {
    coefficients: vec![0; in_out_params.n as usize],
    b: xor,
    current_variance: 0f64,
  };

  temp_result + 2 * ca.clone() + 2 * cb.clone()

  // If the phase is positive, the result is 1/8,
  // otherwise the result is -1/8
  // TODO: Actually implement the bootstrapping so gates can be chained!
  // tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
}

/** bootstrapped Xnor Gate: result = (a==b) */
pub fn boots_xnor(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped Not Gate: result = not(a) */
pub fn boots_not(ca: &LweSample, bk: &TFHEGateBootstrappingCloudKeySet) -> LweSample {
  -ca.clone()
}

/** bootstrapped Copy Gate: result = a */
pub fn boots_copy(ca: LweSample, bk: &TFHEGateBootstrappingCloudKeySet) -> LweSample {
  unimplemented!()
}

/** bootstrapped Nor Gate: result = not(a or b) */
pub fn boots_nor(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  let mu = mod_switch_to_torus32(-1, 8);
  let in_out_params = &bk.params.in_out_params;

  // Compute: (0,-1/8) - ca - cb
  let nor = mod_switch_to_torus32(-1, 8);

  let temp_result = LweSample {
    coefficients: vec![0; in_out_params.n as usize],
    b: nor,
    current_variance: 0f64,
  };

  temp_result - ca.clone() - cb.clone()

  // If the phase is positive, the result is 1/8,
  // otherwise the result is -1/8
  // TODO: Actually implement the bootstrapping so gates can be chained!
  // tfhe_bootstrap_FFT(result, bk->bkFFT, MU, temp_result);
}
/** bootstrapped AndNY Gate: not(a) and b */
pub fn boots_andny(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped AndYN Gate: a and not(b) */
pub fn boots_andyn(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped OrNY Gate: not(a) or b */
pub fn boots_orny(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped OrYN Gate: a or not(b) */
pub fn boots_oryn(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}

/** bootstrapped Mux(a,b,c) = a?b:c */
pub fn boots_mux(
  a: &LweSample,
  b: &LweSample,
  c: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
