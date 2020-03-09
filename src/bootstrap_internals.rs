use crate::lwe::{lwe_key_switch, LweBootstrappingKey, LweSample};
use crate::numerics::{mod_switch_to_torus32, torus_polynomial_mul_by_xai};
use crate::tgsw::{tgsw_extern_mul_to_tlwe, TGswParams, TGswSample};
use crate::tlwe::{tlwe_add_to, tlwe_mul_by_xai_minus_one, TLweSample, Torus32, TorusPolynomial};

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
  bk: &Vec<TGswSample>,
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
  let acc = tfhe_blind_rotate(acc, bk, bara, n, bk_params);
  acc.extract_lwe(extract_params, accum_params)
}

/// Multiply the accumulator by X^sum(bara_i.s_i)
/// # Arguments
/// * `accum` - The TLWE sample to multiply
/// * `bk` - An array of n TGSW samples where bk_i encodes s_i
/// * `bara` - An array of n coefficients between 0 and 2N-1
/// * `bk_params` - The parameters of bk
pub(crate) fn tfhe_blind_rotate(
  accum: TLweSample,
  bk: &Vec<TGswSample>,
  bara: Vec<i32>,
  n: i32,
  bk_params: &TGswParams,
) -> TLweSample {
  let mut temp = TLweSample::new(&bk_params.tlwe_params);
  let temp2 = &mut temp;
  let mut temp3 = accum;

  for i in 0..n as usize {
    let barai = bara[i];
    if barai == 0 {
      // Indeed, this is an easy case!
      // TODO: Figure out why this is an easy case. I'm guessing no work needs to be done for some reason
      continue;
    }
    tfhe_mux_rotate(temp2, &temp3, &bk[i], barai, bk_params);
    std::mem::swap(temp2, &mut temp3);
  }

  // TOOD: Figure out if this was a correct translation
  temp3
}

fn tfhe_mux_rotate(
  result: &mut TLweSample,
  accum: &TLweSample,
  bki: &TGswSample,
  barai: i32,
  bk_params: &TGswParams,
) {
  // ACC = BKi*[(X^barai-1)*ACC]+ACC
  // temp = (X^barai-1)*ACC
  tlwe_mul_by_xai_minus_one(result, barai, accum, &bk_params.tlwe_params);
  // temp *= BKi
  tgsw_extern_mul_to_tlwe(result, bki, bk_params);
  // ACC += temp
  tlwe_add_to(result, accum);
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::tlwe::TLweParameters;
  use rand::Rng;
  const N: u32 = 500;
  const ALPHA_IN: f64 = 5e-4;
  const ALPHA_BK: f64 = 9e-9;

  fn generate_random_key(n: u32) -> Vec<bool> {
    let mut rng = rand::thread_rng();
    (0..n).map(|_| rng.gen()).collect()
  }

  #[test]
  #[ignore]
  fn test_blind_rotate() {
    let mut rng = rand::thread_rng();
    let key = generate_random_key(N);
    let bara: Vec<i32> = (0..N)
      .map(|_| (rng.gen::<u32>() % (2 * N)) as i32)
      .collect();

    let accum_params = TLweParameters::new(N as i32, 1, ALPHA_BK, 1f64 / 16f64);

    let expected_accum_message = TorusPolynomial::torus_polynomial_uniform(N as i32);
    let mut init_alpha_accum = 0.2;
    unimplemented!()
  }
}
