use crate::lwe::{LweBootstrappingKeyFFT, LweSample};
use crate::numerics::mod_switch_to_torus32;
use crate::tlwe::{Torus32, TorusPolynomial};

/// # Arguments
/// * `bk` - The bootstrapping + keyswitch key
/// * `mu` - The output message (if phase(x)>0)
/// * `x` - The input sample
/// returns = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
///
pub(crate) fn bootstrap_fft(bk: &LweBootstrappingKeyFFT, mu: Torus32, x: &LweSample) -> LweSample {
  // LweSample *u = new_LweSample(&bk->accum_params->extracted_lweparams);

  // tfhe_bootstrap_woKS_FFT(u, bk, mu, x);
  // Key switching
  // lweKeySwitch(result, bk->ks, u);

  // delete_LweSample(u);
  unimplemented!()
}

/// # Arguments
/// * `bk` - The bootstrapping + keyswitch key
/// * `mu` - The output message (if phase(x)>0)
/// * `x` - The input sample
/// returns = LWE(mu) iff phase(x)>0, LWE(-mu) iff phase(x)<0
pub(crate) fn bootstrap_without_key_switching(
  bk: &LweBootstrappingKeyFFT,
  mu: Torus32,
  x: &LweSample,
) -> LweSample {
  let bk_params = &bk.bk_params;
  let accum_params = &bk.accum_params;
  let in_params = &bk.in_out_params;
  let big_n = accum_params.n;
  let n_x_2 = 2 * big_n;
  let n = in_params.n;
  let test_vec = TorusPolynomial::new(big_n);
  let mut bara = vec![0; big_n as usize];

  // Modulus switching
  let barb = mod_switch_to_torus32(x.b, n_x_2);
  for i in 0..n as usize {
    bara[i] = mod_switch_to_torus32(x.coefficients[i], n_x_2);
  }

  // The initial testvec = [mu,mu,mu,...,mu]
  // for (int32_t i = 0; i < N; i++) testvect->coefsT[i] = mu;

  // Bootstrapping rotation and extraction
  // tfhe_blindRotateAndExtract_FFT(result, testvect, bk->bkFFT, barb, bara, n, bk_params);

  // delete[] bara;
  // delete_TorusPolynomial(testvect);
  unimplemented!()
}
