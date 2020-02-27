use crate::lwe::{
  LweParams, LweSample, TFHEGateBootstrappingCloudKeySet, TFHEGateBootstrappingParameterSet,
  TFheGateBootstrappingSecretKeySet,
};
use crate::tlwe::TLweParameters;
use crate::tsgw::TGswParams;

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
  //let lwe_key: LweKey = LweKey::new(params.in_out_params.clone());
  // lweKeyGen(lwe_key);
  // TGswKey *tgsw_key = new_TGswKey(params->tgsw_params);
  // tGswKeyGen(tgsw_key);
  // LweBootstrappingKey *bk = new_LweBootstrappingKey(params->ks_t, params->ks_basebit, params->in_out_params,
  //                                                   params->tgsw_params);
  // tfhe_createLweBootstrappingKey(bk, lwe_key, tgsw_key);
  // LweBootstrappingKeyFFT *bkFFT = new_LweBootstrappingKeyFFT(bk);
  // return new TFheGateBootstrappingSecretKeySet(params, bk, bkFFT, lwe_key, tgsw_key);
  unimplemented!()
}

//  deletes gate bootstrapping parameters
// pub fn void delete_gate_bootstrapping_parameters(TFHEGateBootstrappingParameterSet *params);

//  deletes a gate bootstrapping secret key
//pub fn void delete_gate_bootstrapping_secret_keyset(TFheGateBootstrappingSecretKeySet *keyset);

//  deletes a gate bootstrapping secret key
//pub fn void delete_gate_bootstrapping_cloud_keyset(TFHEGateBootstrappingCloudKeySet *keyset);

/** generate a new unititialized ciphertext */
pub fn new_gate_bootstrapping_ciphertext(params: &TFHEGateBootstrappingParameterSet) -> LweSample {
  unimplemented!()
}

/** generate a new unititialized ciphertext array of length nbelems */
pub fn new_gate_bootstrapping_ciphertext_array(
  nbelems: i32,
  params: &TFHEGateBootstrappingParameterSet,
) -> LweSample {
  unimplemented!()
}

//  deletes a ciphertext
// pub fn void delete_gate_bootstrapping_ciphertext(LweSample *sample);

//  deletes a ciphertext array of length nbelems
// pub fn void delete_gate_bootstrapping_ciphertext_array(int32_t nbelems, LweSample *samples);

/** encrypts a boolean */
pub fn boots_sym_encrypt(message: bool, params: &TFheGateBootstrappingSecretKeySet) -> LweSample {
  unimplemented!()
}

/** decrypts a boolean */
pub fn boots_sym_decrypt(sample: &LweSample, params: &TFheGateBootstrappingSecretKeySet) -> bool {
  unimplemented!()
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
