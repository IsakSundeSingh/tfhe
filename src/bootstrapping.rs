use crate::lwe::{
  LweKey, LweParams, LweSample, TFHEGateBootstrappingCloudKeySet,
  TFHEGateBootstrappingParameterSet, TFheGateBootstrappingSecretKeySet,
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
  const k: i32 = 1;
  const n: i32 = 500;
  const bk_l: i32 = 2;
  const bk_Bgbit: i32 = 10;
  const ks_basebit: i32 = 2;
  const ks_length: i32 = 8;
  const ks_stdev: f64 = 2.44e-5; //standard deviation
  const bk_stdev: f64 = 7.18e-9; //standard deviation
  const max_stdev: f64 = 0.012467; //max standard deviation for a 1/4 msg space
  let params_in: LweParams = LweParams::new(n, ks_stdev, max_stdev);
  let params_accum: TLweParameters = TLweParameters::new(N, k, bk_stdev, max_stdev);
  let params_bk: TGswParams = TGswParams::new(bk_l, bk_Bgbit, params_accum);

  // TfheGarbageCollector::register_param(params_in);
  // TfheGarbageCollector::register_param(params_accum);
  // TfheGarbageCollector::register_param(params_bk);
  TFHEGateBootstrappingParameterSet::new(ks_length, ks_basebit, params_in, params_bk)
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
pub fn bootsSymEncrypt(message: bool, params: &TFheGateBootstrappingSecretKeySet) -> LweSample {
  unimplemented!()
}

/** decrypts a boolean */
pub fn bootsSymDecrypt(sample: &LweSample, params: &TFheGateBootstrappingSecretKeySet) -> bool {
  unimplemented!()
}

/** bootstrapped Constant (true or false) trivial Gate */
pub fn bootsCONSTANT(value: i32, bk: &TFHEGateBootstrappingCloudKeySet) -> LweSample {
  unimplemented!()
}

/** bootstrapped Nand Gate */
pub fn bootsNAND(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped Or Gate:  */
pub fn bootsOR(ca: &LweSample, cb: &LweSample, bk: &TFHEGateBootstrappingCloudKeySet) -> LweSample {
  unimplemented!()
}
/** bootstrapped And Gate: result = a and b */
pub fn bootsAND(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped Xor Gate: result = a xor b */
pub fn bootsXOR(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped Xnor Gate: result = (a==b) */
pub fn bootsXNOR(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped Not Gate: result = not(a) */
pub fn bootsNOT(ca: &LweSample, bk: &TFHEGateBootstrappingCloudKeySet) -> LweSample {
  unimplemented!()
}

/** bootstrapped Copy Gate: result = a */
pub fn bootsCOPY(ca: LweSample, bk: &TFHEGateBootstrappingCloudKeySet) -> LweSample {
  unimplemented!()
}
/** bootstrapped Nor Gate: result = not(a or b) */
pub fn bootsNOR(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped AndNY Gate: not(a) and b */
pub fn bootsANDNY(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped AndYN Gate: a and not(b) */
pub fn bootsANDYN(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped OrNY Gate: not(a) or b */
pub fn bootsORNY(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
/** bootstrapped OrYN Gate: a or not(b) */
pub fn bootsORYN(
  ca: &LweSample,
  cb: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}

/** bootstrapped Mux(a,b,c) = a?b:c */
pub fn bootsMUX(
  a: &LweSample,
  b: &LweSample,
  c: &LweSample,
  bk: &TFHEGateBootstrappingCloudKeySet,
) -> LweSample {
  unimplemented!()
}
