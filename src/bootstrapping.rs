use crate::lwe::LweKey;
use crate::lwe::{
  LweBootstrappingKey, LweParams, LweSample, TFHEGateBootstrappingCloudKeySet,
  TFHEGateBootstrappingParameterSet, TFheGateBootstrappingSecretKeySet,
};
use crate::tgsw::{TGswKey, TGswParams};
use crate::tlwe::TLweParameters;

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
  const KS_STDEV: f64 = 2.44e-5; //standard deviation
  const BK_STDEV: f64 = 7.18e-9; //standard deviation
  const MAX_STDEV: f64 = 0.012_467; //max standard deviation for a 1/4 msg space
  let params_in: LweParams = LweParams::new(LOWERCASE_N, KS_STDEV, MAX_STDEV);
  let params_accum: TLweParameters = TLweParameters::new(N, K, BK_STDEV, MAX_STDEV);
  let params_bk: TGswParams = TGswParams::new(BK_L, BK_BG_BIT, params_accum);

  // TfheGarbageCollector::register_param(params_in);
  // TfheGarbageCollector::register_param(params_accum);
  // TfheGarbageCollector::register_param(params_bk);
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
  // LweBootstrappingKeyFFT *bkFFT = new_LweBootstrappingKeyFFT(bk);
  // return new TFheGateBootstrappingSecretKeySet(params, bk, bkFFT, lwe_key, tgsw_key);
}

//  deletes gate bootstrapping parameters
// pub fn void delete_gate_bootstrapping_parameters(TFHEGateBootstrappingParameterSet *params);

//  deletes a gate bootstrapping secret key
//pub fn void delete_gate_bootstrapping_secret_keyset(TFheGateBootstrappingSecretKeySet *keyset);

//  deletes a gate bootstrapping secret key
//pub fn void delete_gate_bootstrapping_cloud_keyset(TFHEGateBootstrappingCloudKeySet *keyset);

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

//  deletes a ciphertext
// pub fn void delete_gate_bootstrapping_ciphertext(LweSample *sample);

//  deletes a ciphertext array of length nbelems
// pub fn void delete_gate_bootstrapping_ciphertext_array(int32_t nbelems, LweSample *samples);

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
  unimplemented!()
}
/** bootstrapped Or Gate:  */
pub fn boots_or(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped And Gate: result = a and b */
pub fn boots_and(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped Xor Gate: result = a xor b */
pub fn boots_xor(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
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
  unimplemented!()
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
  unimplemented!()
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
