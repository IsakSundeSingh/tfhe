//! Module for all homomorphic binary and unary gate operations.

#[cfg(feature = "bootstrapping")]
use crate::bootstrapping::tfhe_bootstrap;

use crate::{
  lwe::{CloudKey, LweSample, Parameters},
  numerics::{encode_message, Torus32},
};

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
  let in_out_params = &bk.params.in_out_params;

  // Compute: (0,1/8) - ca - cb
  const NAND: Torus32 = encode_message(1, 8);
  let temp_result = LweSample::trivial(NAND, in_out_params);

  // If the phase is positive, the result is 1/8,
  // otherwise the result is -1/8
  #[cfg(feature = "bootstrapping")]
  {
    const MU: Torus32 = encode_message(1, 8);
    tfhe_bootstrap(&bk.bk, MU, temp_result - ca - cb)
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result - ca - cb
  }
}

/** bootstrapped Or Gate:  */
pub fn boots_or(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  let in_out_params = &bk.params.in_out_params;

  // Compute: (0,1/8) + ca + cb
  const OR: Torus32 = encode_message(1, 8);

  let temp_result = LweSample::trivial(OR, in_out_params);

  #[cfg(feature = "bootstrapping")]
  {
    const MU: Torus32 = encode_message(1, 8);
    tfhe_bootstrap(&bk.bk, MU, temp_result + ca + cb)
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result + ca + cb
  }
}

/** bootstrapped And Gate: result = a and b */
pub fn boots_and(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  let in_out_params = &bk.params.in_out_params;

  const AND: Torus32 = encode_message(-1, 8);

  // Compute: (0,-1/8) + ca + cb
  let temp_result = LweSample::trivial(AND, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    const MU: Torus32 = encode_message(1, 8);
    tfhe_bootstrap(&bk.bk, MU, temp_result + ca + cb)
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result + ca + cb
  }
}

/** bootstrapped Xor Gate: result = a xor b */
pub fn boots_xor(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  let in_out_params = &bk.params.in_out_params;

  const XOR: Torus32 = encode_message(1, 4);

  // Compute: (0,1/4) + 2*(ca + cb)
  let temp_result = LweSample::trivial(XOR, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    const MU: Torus32 = encode_message(1, 8);
    tfhe_bootstrap(&bk.bk, MU, temp_result + 2 * ca + 2 * cb)
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result + 2 * ca + 2 * cb
  }
}

/** bootstrapped Xnor Gate: result = (a==b) */
pub fn boots_xnor(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  // use crate::bootstrap_internals::tfhe_bootstrap;
  let in_out_params = &bk.params.in_out_params;

  const XNOR: Torus32 = encode_message(-1, 4);

  // Compute: (0,-1/4) + 2*(-ca-cb)
  let temp_result = LweSample::trivial(XNOR, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    const MU: Torus32 = encode_message(1, 8);
    tfhe_bootstrap(&bk.bk, MU, temp_result - 2 * ca - 2 * cb)
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result - 2 * ca - 2 * cb
  }
}

/** bootstrapped Not Gate: result = not(a) */
pub fn boots_not(ca: &LweSample, _bk: &CloudKey) -> LweSample {
  !ca
}

/** bootstrapped Nor Gate: result = not(a or b) */
pub fn boots_nor(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  let in_out_params = &bk.params.in_out_params;

  const NOR: Torus32 = encode_message(-1, 8);

  // Compute: (0,-1/8) - ca - cb
  let temp_result = LweSample::trivial(NOR, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    const MU: Torus32 = encode_message(-1, 8);
    tfhe_bootstrap(&bk.bk, MU, temp_result - ca - cb)
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result - ca - cb
  }
}

/** bootstrapped AndNY Gate: not(a) and b */
pub fn boots_andny(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  let in_out_params = &bk.params.in_out_params;

  // Compute: (0,-1/8) - ca + cb
  const ANDNY: Torus32 = encode_message(-1, 8);
  let temp_result = LweSample::trivial(ANDNY, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    const MU: Torus32 = encode_message(1, 8);
    tfhe_bootstrap(&bk.bk, MU, temp_result - ca + cb)
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result - ca + cb
  }
}

/** bootstrapped AndYN Gate: a and not(b) */
pub fn boots_andyn(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  let in_out_params = &bk.params.in_out_params;

  const ANDYN: Torus32 = encode_message(-1, 8);

  // Compute: (0,-1/8) + ca - cb
  let temp_result = LweSample::trivial(ANDYN, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    const MU: Torus32 = encode_message(1, 8);
    tfhe_bootstrap(&bk.bk, MU, temp_result + ca - cb)
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result + ca - cb
  }
}

/** bootstrapped OrNY Gate: not(a) or b */
pub fn boots_orny(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  let in_out_params = &bk.params.in_out_params;

  const ORNY: Torus32 = encode_message(1, 8);

  // Compute: (0,1/8) - ca + cb
  let temp_result = LweSample::trivial(ORNY, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    const MU: Torus32 = encode_message(1, 8);
    tfhe_bootstrap(&bk.bk, MU, temp_result - ca + cb)
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result - ca + cb
  }
}

/** bootstrapped OrYN Gate: a or not(b) */
pub fn boots_oryn(ca: &LweSample, cb: &LweSample, bk: &CloudKey) -> LweSample {
  let in_out_params = &bk.params.in_out_params;

  const ORYN: Torus32 = encode_message(1, 8);

  // Compute: (0,1/8) + ca - cb
  let temp_result = LweSample::trivial(ORYN, in_out_params);
  #[cfg(feature = "bootstrapping")]
  {
    const MU: Torus32 = encode_message(1, 8);
    tfhe_bootstrap(&bk.bk, MU, temp_result + ca - cb)
  }
  #[cfg(not(feature = "bootstrapping"))]
  {
    temp_result + ca - cb
  }
}

/// Bootstrapped MUX-gate.
/// Essentially performs the following
/// ```rust,no_run
/// # let a = true;
/// # let b = true;
/// # let c = true;
/// # let _ =
/// if a {
///   b
/// } else {
///   c
/// }
/// # ;
/// ```
#[allow(unused)]
pub fn boots_mux(a: &LweSample, b: &LweSample, c: &LweSample, bk: &CloudKey) -> LweSample {
  #[cfg(not(feature = "bootstrapping"))]
  {
    panic!("Mux gate cannot run without bootstrapping");
  }

  #[cfg(feature = "bootstrapping")]
  {
    use crate::bootstrapping::tfhe_bootstrap_without_key_switching;
    use crate::lwe::lwe_key_switch;

    const MU: Torus32 = encode_message(1, 8);
    let in_out_params = &bk.params.in_out_params;
    let extracted_params = &bk.params.tgsw_params.tlwe_params.extracted_lweparams;

    const AND: Torus32 = encode_message(-1, 8);

    // Compute AND(a,b) = (0, -1/8) + a + b
    let temp_result = LweSample::trivial(AND, in_out_params) + a + b;
    let u1 = tfhe_bootstrap_without_key_switching(&bk.bk, MU, temp_result);

    // Compute AND(NOT(a), c) = (0, -1/8) - a + c
    let temp_result = LweSample::trivial(AND, in_out_params) - a + c;
    let u2 = tfhe_bootstrap_without_key_switching(&bk.bk, MU, temp_result);

    const MUX: Torus32 = encode_message(1, 8);
    let temp_result = LweSample::trivial(MUX, extracted_params) + u1 + u2;
    lwe_key_switch(&bk.bk.ks, temp_result)
  }
}
