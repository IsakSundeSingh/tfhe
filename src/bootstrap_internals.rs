use crate::lwe::{LweBootstrappingKey, LweSample};
use crate::numerics::mod_switch_to_torus32;
use crate::{
  tgsw::{TGswParams, TGswSample},
  tlwe::{TLweSample, Torus32, TorusPolynomial},
};

/// # Arguments
/// * `bk` - The bootstrapping + keyswitch key
/// * `mu` - The output message (if phase(x)>0)
/// * `x` - The input sample
/// returns = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
pub(crate) fn tfhe_bootstrap(bk: &LweBootstrappingKey, mu: Torus32, x: &LweSample) -> LweSample {
  // LweSample *u = new_LweSample(&bk->accum_params->extracted_lweparams);

  // tfhe_bootstrap_woKS(u, bk, mu, x);
  // Key Switching
  // lweKeySwitch(result, bk->ks, u);

  // delete_LweSample(u);
  unimplemented!()
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
  x: &LweSample,
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

  // tfhe_blindRotateAndExtract(result, testvect, bk->bk, barb, bara, n, bk_params);

  // delete[] bara;
  // delete_TorusPolynomial(testvect);
  unimplemented!()
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
  bk: &TGswSample,
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
  // tfhe_blindRotate(acc, bk, bara, n, bk_params);
  // tLweExtractLweSample(result, acc, extract_params, accum_params);

  // delete_TLweSample(acc);
  // delete_TorusPolynomial(testvectbis);
  unimplemented!()
}

//result= X^{a}*source
pub(crate) fn torus_polynomial_mul_by_xai(a: i32, source: &TorusPolynomial) -> TorusPolynomial {
  let n = source.n;
  let r#in = &source.coefs;
  assert!(a >= 0 && a < 2 * n);
  let mut coefs = vec![0; n as usize];

  if a < n {
    for i in 0..a {
      // So that i-a<0 (French: sur que ...)
      coefs[i as usize] = -source.coefs[(i - a + n) as usize];
    }
    for i in a..n {
      // So that N>i-a>=0 (French: sur que ...)
      coefs[i as usize] = source.coefs[(i - a) as usize];
    }
  } else {
    let aa = a - n;
    for i in 0..aa {
      // So that i-a<0 (French: sur que ...)
      coefs[i as usize] = source.coefs[(i - aa + n) as usize];
    }
    for i in aa..n {
      // So that N>i-a>=0 (French: sur que ...)
      coefs[i as usize] = -source.coefs[(i - aa) as usize];
    }
  }

  let n = coefs.len() as i32;
  TorusPolynomial { coefs, n }
}

/**
 * multiply the accumulator by X^sum(bara_i.s_i)
 * @param accum the TLWE sample to multiply
 * @param bk An array of n TGSW samples where bk_i encodes s_i
 * @param bara An array of n coefficients between 0 and 2N-1
 * @param bk_params The parameters of bk
 */
pub(crate) fn tfhe_blind_rotate(
  accum: &mut TLweSample,
  bk: Vec<TGswSample>,
  bara: Vec<i32>,
  n: i32,
  bk_params: &TGswParams,
) {
  let temp = TLweSample::new(&bk_params.tlwe_params);
  let temp2 = &temp;
  let temp3 = &accum;
  for i in 0..n as usize {
    let barai = bara[i];
    if barai == 0 {
      // Indeed, this is an easy case!
      continue;
    }

    //tfhe_MuxRotate(temp2, temp3, bk + i, barai, bk_params);
    //swap(temp2, temp3);
  }
  // if (temp3 != accum) {
  //     tLweCopy(accum, temp3, bk_params->tlwe_params);
  // }

  // delete_TLweSample(temp);

  unimplemented!()
}
