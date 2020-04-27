use crate::{
  lwe::{CloudKey, LweBootstrappingKey, LweKey, LweSample, Parameters, SecretKey},
  numerics::{encode_message, Torus32},
  tgsw::TGswKey,
};

#[cfg(feature = "bootstrapping")]
use crate::bootstrap_internals::tfhe_bootstrap;

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
/// # use tfhe::bootstrapping::{generate_keys, bootstrapping_parameters};
/// # use tfhe::lwe::Parameters;
/// let params = Parameters::default();
/// let (secret_key, cloud_key) = generate_keys(&params);
/// ```
pub fn generate_keys(params: &Parameters) -> (SecretKey, CloudKey) {
  let lwe_key: LweKey = LweKey::generate(&params.in_out_params);
  let tgsw_key = TGswKey::generate(&params.tgsw_params);
  let bk = LweBootstrappingKey::create(&params, &lwe_key, &tgsw_key);
  (
    SecretKey::new(params.clone(), lwe_key, tgsw_key),
    CloudKey::new(params.clone(), bk),
  )
}

/** generate a new unititialized ciphertext */
pub fn new_gate_bootstrapping_ciphertext(params: &Parameters) -> LweSample {
  LweSample::new(&params.in_out_params)
}

/** generate a new unititialized ciphertext array of length nbelems */
pub fn new_gate_bootstrapping_ciphertext_array(
  nbelems: i32,
  params: &Parameters,
) -> Vec<LweSample> {
  vec![new_gate_bootstrapping_ciphertext(&params); nbelems as usize]
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

/// Convert some constant to an encoded form so it can be used in the homomorphic operations with the ciphertexts.///
pub fn boots_constant(value: bool, bk: &CloudKey) -> LweSample {
  let in_out_params = &bk.params.in_out_params;
  const MU: Torus32 = encode_message(1, 8);
  LweSample {
    coefficients: vec![0; in_out_params.n as usize],
    b: if value { MU } else { -MU },
    current_variance: 0_f64,
  }
}

/** bootstrapped Nand Gate */
pub fn boots_nand(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  const MU: Torus32 = encode_message(1, 8);
  let in_out_params = &bk.params.in_out_params;

  // Compute: (0,1/8) - ca - cb
  const NAND: Torus32 = encode_message(1, 8);
  let temp_result = LweSample::trivial(NAND, in_out_params);

  // If the phase is positive, the result is 1/8,
  // otherwise the result is -1/8
  // TODO: Actually implement the bootstrapping so gates can be chained!
  #[cfg(feature = "bootstrapping")]
  {
    tfhe_bootstrap(&bk.bk, MU, temp_result - ca.clone() - cb.clone())
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result - ca.clone() - cb.clone()
  }
}

/** bootstrapped Or Gate:  */
pub fn boots_or(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  const MU: Torus32 = encode_message(1, 8);
  let in_out_params = &bk.params.in_out_params;

  // Compute: (0,1/8) + ca + cb
  const OR: Torus32 = encode_message(1, 8);

  let temp_result = LweSample::trivial(OR, in_out_params);

  #[cfg(feature = "bootstrapping")]
  {
    tfhe_bootstrap(&bk.bk, MU, temp_result + ca.clone() + cb.clone())
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result + ca.clone() + cb.clone()
  }
}

/** bootstrapped And Gate: result = a and b */
pub fn boots_and(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  const MU: Torus32 = encode_message(1, 8);
  let in_out_params = &bk.params.in_out_params;

  const AND: Torus32 = encode_message(-1, 8);

  // Compute: (0,-1/8) + ca + cb
  let temp_result = LweSample::trivial(AND, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    tfhe_bootstrap(&bk.bk, MU, temp_result + ca.clone() + cb.clone())
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result + ca.clone() + cb.clone()
  }
}

/** bootstrapped Xor Gate: result = a xor b */
pub fn boots_xor(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  const MU: Torus32 = encode_message(1, 8);
  let in_out_params = &bk.params.in_out_params;

  const XOR: Torus32 = encode_message(1, 4);

  // Compute: (0,1/4) + 2*(ca + cb)
  let temp_result = LweSample::trivial(XOR, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    tfhe_bootstrap(&bk.bk, MU, temp_result + 2 * ca.clone() + 2 * cb.clone())
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result + 2 * ca.clone() + 2 * cb.clone()
  }
}

/** bootstrapped Xnor Gate: result = (a==b) */
pub fn boots_xnor(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  // use crate::bootstrap_internals::tfhe_bootstrap;
  const MU: Torus32 = encode_message(1, 8);
  let in_out_params = &bk.params.in_out_params;

  const XNOR: Torus32 = encode_message(-1, 4);

  // Compute: (0,-1/4) + 2*(-ca-cb)
  let temp_result = LweSample::trivial(XNOR, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    tfhe_bootstrap(&bk.bk, MU, temp_result - 2 * ca.clone() - 2 * cb.clone())
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result - 2 * ca.clone() - 2 * cb.clone()
  }
}

/** bootstrapped Not Gate: result = not(a) */
pub fn boots_not(ca: &LweSample, _bk: &CloudKey) -> LweSample {
  !ca.clone()
}

/** bootstrapped Nor Gate: result = not(a or b) */
pub fn boots_nor(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  const MU: Torus32 = encode_message(-1, 8);
  let in_out_params = &bk.params.in_out_params;

  const NOR: Torus32 = encode_message(-1, 8);

  // Compute: (0,-1/8) - ca - cb
  let temp_result = LweSample::trivial(NOR, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    tfhe_bootstrap(&bk.bk, MU, temp_result - ca.clone() - cb.clone())
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result - ca.clone() - cb.clone()
  }
}

/** bootstrapped AndNY Gate: not(a) and b */
pub fn boots_andny(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  const MU: Torus32 = encode_message(1, 8);
  let in_out_params = &bk.params.in_out_params;

  // Compute: (0,-1/8) - ca + cb
  const ANDNY: Torus32 = encode_message(-1, 8);
  let temp_result = LweSample::trivial(ANDNY, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    tfhe_bootstrap(&bk.bk, MU, temp_result - ca.clone() + cb.clone())
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result - ca.clone() + cb.clone()
  }
}

/** bootstrapped AndYN Gate: a and not(b) */
pub fn boots_andyn(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  const MU: Torus32 = encode_message(1, 8);
  let in_out_params = &bk.params.in_out_params;

  const ANDYN: Torus32 = encode_message(-1, 8);

  // Compute: (0,-1/8) + ca - cb
  let temp_result = LweSample::trivial(ANDYN, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    tfhe_bootstrap(&bk.bk, MU, temp_result + ca.clone() - cb.clone())
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result + ca.clone() - cb.clone()
  }
}

/** bootstrapped OrNY Gate: not(a) or b */
pub fn boots_orny(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  const MU: Torus32 = encode_message(1, 8);
  let in_out_params = &bk.params.in_out_params;

  const ORNY: Torus32 = encode_message(1, 8);

  // Compute: (0,1/8) - ca + cb
  let temp_result = LweSample::trivial(ORNY, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    tfhe_bootstrap(&bk.bk, MU, temp_result - ca.clone() + cb.clone())
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result - ca.clone() + cb.clone()
  }
}

/** bootstrapped OrYN Gate: a or not(b) */
pub fn boots_oryn(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  const MU: Torus32 = encode_message(1, 8);
  let in_out_params = &bk.params.in_out_params;

  const ORYN: Torus32 = encode_message(1, 8);

  // Compute: (0,1/8) + ca - cb
  let temp_result = LweSample::trivial(ORYN, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    tfhe_bootstrap(&bk.bk, MU, temp_result + ca.clone() - cb.clone())
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result + ca.clone() - cb.clone()
  }
}

/** bootstrapped Mux(a,b,c) = a?b:c */
pub fn boots_mux(_a: &LweSample, _b: &LweSample, _c: &LweSample, _bk: &CloudKey) -> LweSample {
  todo!()
}
