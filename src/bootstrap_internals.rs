use crate::lwe::{lwe_key_switch, LweBootstrappingKey, LweSample};
use crate::numerics::{mod_switch_to_torus32, torus_polynomial_mul_by_xai, Torus32};
use crate::polynomial::TorusPolynomial;
use crate::tgsw::{tgsw_extern_mul_to_tlwe, TGswParams, TGswSample};
use crate::tlwe::{tlwe_add_to, tlwe_mul_by_xai_minus_one, TLweSample};

/// # Arguments
/// * `bk` - The bootstrapping + keyswitch key
/// * `mu` - The output message (if phase(x)>0)
/// * `x` - The input sample
/// returns = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
pub(crate) fn tfhe_bootstrap(bk: &LweBootstrappingKey, mu: Torus32, x: LweSample) -> LweSample {
  let res = tfhe_bootstrap_without_key_switching(bk, mu, x);
  // Key Switching
  lwe_key_switch(&bk.ks, res)
}

/**
 * result = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
 * @param result The resulting LweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
 */
pub(crate) fn tfhe_bootstrap_without_key_switching(
  bk: &LweBootstrappingKey,
  mu: Torus32,
  x: LweSample,
) -> LweSample {
  let bk_params = &bk.bk_params;
  let accum_params = &bk.accum_params;
  let in_params = &bk.in_out_params;
  let big_n = accum_params.n;
  let n_x_2 = 2 * big_n;
  let n = in_params.n;
  let mut test_vec = TorusPolynomial::new(big_n);
  let mut bara = vec![0; big_n as usize];

  let barb = mod_switch_to_torus32(x.b, n_x_2);
  for i in 0..n as usize {
    bara[i] = mod_switch_to_torus32(x.coefficients[i], n_x_2);
  }

  // The initial testvec = [mu,mu,mu,...,mu]
  for i in 0..big_n as usize {
    test_vec.coefs[i] = mu;
  }

  tfhe_blind_rotate_and_extract(test_vec, &bk.bk, barb, bara, n, bk_params)
}

/**
 * result = LWE(v_p) where p=barb-sum(bara_i.s_i) mod 2N
 * @param result the output LWE sample
 * @param v a 2N-elt anticyclic function (represented by a TorusPolynomial)
 * @param bk An array of n TGSW samples where bk_i encodes s_i
 * @param barb A coefficients between 0 and 2N-1
 * @param bara An array of n coefficients between 0 and 2N-1
 * @param bk_params The parameters of bk
 */

pub(crate) fn tfhe_blind_rotate_and_extract(
  v: TorusPolynomial,
  bk: &[TGswSample],
  barb: i32,
  bara: Vec<i32>,
  n: i32,
  bk_params: &TGswParams,
) -> LweSample {
  let accum_params = &bk_params.tlwe_params;
  let extract_params = &accum_params.extracted_lweparams;
  let n = accum_params.n;
  let _2n = 2 * n;
  let test_vec_bis = if barb != 0 {
    torus_polynomial_mul_by_xai(_2n - barb, &v)
  } else {
    v
  };

  // if (barb != 0) torusPolynomialMulByXai(testvectbis, _2N - barb, v);
  // else torusPolynomialCopy(testvectbis, v);

  let acc = TLweSample::trivial(test_vec_bis, accum_params);
  let acc = tfhe_blind_rotate(acc, bk, &bara[..], n, bk_params);
  acc.extract_lwe(extract_params, accum_params)
}

/// The BlindRotate algorithm multiplies the polynomial encrypted in the input
/// TRLWE ciphertext by an encrypted power of X. The effect produced is a rotation of the coefficients.
/// Multiply the accumulator by X^sum(bara_i.s_i)
/// # Arguments
/// * `accum` - The TLWE sample to multiply
/// * `bk` - An array of n TGSW samples where bk_i encodes s_i
/// * `bara` - An array of n coefficients between 0 and 2N-1
/// * `bk_params` - The parameters of bk
pub(crate) fn tfhe_blind_rotate(
  accum: TLweSample,
  bk: &[TGswSample],
  bara: &[i32],
  n: i32,
  bk_params: &TGswParams,
) -> TLweSample {
  let mut temp = TLweSample::new(&bk_params.tlwe_params);
  let mut temp2 = accum;

  for i in 0..n as usize {
    let barai = bara[i];
    if barai == 0 {
      // Indeed, this is an easy case!
      // TODO: Figure out why this is an easy case. I'm guessing no work needs to be done for some reason
      continue;
    }

    temp = tfhe_mux_rotate(&temp.clone(), &temp2, &bk[i], barai, bk_params);
    std::mem::swap(&mut temp, &mut temp2);
  }

  // TODO: Figure out if this was a correct translation
  temp2
}

fn tfhe_mux_rotate(
  result: &TLweSample,
  accum: &TLweSample,
  bki: &TGswSample,
  barai: i32,
  bk_params: &TGswParams,
) -> TLweSample {
  // ACC = BKi*[(X^barai-1)*ACC]+ACC
  // temp = (X^barai-1)*ACC
  let res = tlwe_mul_by_xai_minus_one(result, barai, accum, &bk_params.tlwe_params);
  // temp *= BKi
  let res = tgsw_extern_mul_to_tlwe(&res, bki, bk_params);
  res + accum.clone()
  // ACC += temp
  // tlwe_add_to(result, accum);
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::{polynomial::Polynomial, tlwe::TLweParameters};
  use rand::Rng;
  use rand_distr::Distribution;

  const BIG_N: u32 = 1024;
  const SMOL_N: u32 = 500;
  const ALPHA_IN: f64 = 5e-4;
  const ALPHA_BK: f64 = 9e-9;
  const BG_BIT_BK: i32 = 10;
  const K: i32 = 1;
  const L_BK: i32 = 3;

  fn generate_random_key(n: u32) -> Vec<bool> {
    let mut rng = rand::thread_rng();
    (0..n).map(|_| rng.gen()).collect()
  }

  #[test]
  #[ignore]
  fn test_blind_rotate() {
    let d = rand_distr::Uniform::new(0, i32::max_value());
    let mut rng = rand::thread_rng();
    let key = generate_random_key(SMOL_N);
    let bara: Vec<i32> = (0..SMOL_N)
      /* NOTE: This actually said BIG_N,
      but that doesn't make sense as it's greater than SMOL_N
      and will cause a panic in crate::numerics::torus_polynomial_mul_by_xai_minus_one;*/
      .map(|_| (d.sample(&mut rng) % (2 * SMOL_N as i32)))
      .collect();

    let accum_params = TLweParameters::new(SMOL_N as i32, 1, ALPHA_BK, 1f64 / 16f64);
    let bk_params = TGswParams::new(L_BK, BG_BIT_BK, accum_params.clone());
    let bk: Vec<TGswSample> = (0..SMOL_N).map(|_| TGswSample::new(&bk_params)).collect();

    let mut expected_accum_message = TorusPolynomial::uniform(BIG_N as usize);
    let init_accum_message = expected_accum_message.clone();
    let init_alpha_accum = 0.2;
    let mut accum = TLweSample::new(&accum_params);
    accum.current_variance = init_alpha_accum * init_alpha_accum;
    let mut expected_offset = 0;

    // for i in 0..SMOL_N as usize {
    //   accum = tfhe_blind_rotate(accum, &bk[i..], &bara[i..], 1, &bk_params);

    //   if key[i] == true && bara[i] != 0 {
    //     expected_offset = (expected_offset + bara[i]) % ((2 * BIG_N) as i32);
    //     expected_accum_message = torus_polynomial_mul_by_xai(expected_offset, &init_accum_message);
    //   }

    //   expected_accum_message
    //     .coefs
    //     .iter()
    //     .zip(accum.b.coefs.iter())
    //     .for_each(|(a, b)| assert_eq!(a, b))
    // }

    // Now, bootstraprotate: all iterations at once (same offset)
    // torusPolynomialCopy(faccum->message, initAccumMessage);
    // faccum->current_variance = initAlphaAccum * initAlphaAccum;

    accum = tfhe_blind_rotate(accum, &bk, &bara, SMOL_N as i32, &bk_params);
    expected_accum_message
      .coefs
      .iter()
      .zip(accum.b.coefs.iter())
      .for_each(|(a, b)| assert_eq!(a, b));
  }
}
