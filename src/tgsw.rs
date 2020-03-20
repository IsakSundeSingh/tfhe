use crate::tlwe::{
  IntPolynomial, TLweKey, TLweParameters, TLweSample, TLweSampleFFT, Torus32, TorusPolynomial,
};

#[derive(Clone, Debug, PartialEq)]
pub struct TGswParams {
  /// decomp length
  l: i32,

  /// log_2(Bg)
  bg_bit: i32,

  /// decomposition base (must be a power of 2)
  bg: i32,

  /// Bg/2
  half_bg: i32,

  /// Bg-1
  mask_mod: u32,

  /// Params of each row
  pub(crate) tlwe_params: TLweParameters,

  /// number of rows = (k+1)*l
  kpl: i32,

  /// powers of Bgbit
  h: Vec<Torus32>,

  /// offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))
  offset: u32,
}

impl TGswParams {
  pub fn new(l: i32, bg_bit: i32, tlwe_params: TLweParameters) -> Self {
    let bg = (1 << bg_bit) as i32;
    let half_bg = (bg >> 1) as i32;
    let mask_mod = (bg - 1) as u32;
    let mut h = vec![l];
    for i in 0..l {
      let kk = (32 - (i + 1)) * bg_bit;
      // 1/(Bg^(i+1)) as a Torus32
      h.push(1i32.checked_shl(kk as u32).unwrap_or(0));
    }

    // offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))
    let temp: u32 = (0..1).map(|i| 1 << (32 - (i + 1) * bg_bit)).sum();
    let offset = temp * half_bg as u32;

    let kpl = (tlwe_params.k + 1) * l;

    Self {
      l,
      bg_bit,
      bg,
      half_bg,
      mask_mod,
      tlwe_params,
      kpl,
      h,
      offset,
    }
  }
}

pub struct TGswKey {
  /// the parameters
  pub(crate) params: TGswParams,

  /// the tlwe params of each rows
  tlwe_params: TLweParameters,

  /// the key (array of k polynomials)
  key: Vec<IntPolynomial>,

  pub(crate) tlwe_key: TLweKey,
}

impl TGswKey {
  // same key as in TLwe
  fn new(params: &TGswParams) -> Self {
    let tlwe_params = params.tlwe_params.clone();
    let tlwe_key = TLweKey::new(&tlwe_params);
    let key = tlwe_key.key.clone();
    Self {
      params: params.clone(),
      tlwe_params,
      key,
      tlwe_key,
    }
  }

  pub(crate) fn generate(params: &TGswParams) -> Self {
    let mut key = Self::new(params);
    key.tlwe_key.generate();
    key
  }

  pub(crate) fn encrypt(&self, result: &mut TGswSample, message: i32, alpha: f64) {
    result.encrypt_zero(alpha, &self);
    result.add_mu_int_h(message, &self.params);
  }
}

#[derive(Clone, Debug, PartialEq)]
pub struct TGswSample {
  /// (k+1)l TLwe Sample (THIS IS A MATRIX)
  all_sample: Vec<Vec<TLweSample>>,
  /// Horizontal blocks (l lines) of TGSW matrix
  // bloc_sample: Vec<TLweSample>,
  k: i32,
  l: i32,
}

impl TGswSample {
  pub(crate) fn new(params: &TGswParams) -> Self {
    let k = params.tlwe_params.k;
    // Lines / rows
    let l = params.l;
    // TODO: find out if this is correctamente
    let all_sample = vec![vec![TLweSample::new(&params.tlwe_params); (k + 1) as usize]; l as usize];

    Self { all_sample, k, l }
  }

  pub(crate) fn encrypt_zero(&mut self, alpha: f64, key: &TGswKey) {
    let rl_key = &key.tlwe_key;
    let kpl = key.params.kpl;

    self.all_sample[0]
      .iter_mut()
      .for_each(|s| s.encrypt_zero(alpha, rl_key));
    // for p in 0..kpl as usize {
    //   self.all_sample[0][p].encrypt_zero(alpha, rl_key);
    // }
  }

  pub(crate) fn add_mu_int_h(&mut self, message: i32, params: &TGswParams) {
    let k = params.tlwe_params.k;
    let l = params.l;
    let h = &params.h;

    // TFHE comment: Compute self += H
    // My comment:   Compute self += H * message (ish)
    let hs: Vec<i32> = h.iter().map(|x| message * x).collect();
    self.all_sample = self
      .all_sample
      .iter()
      .enumerate()
      .map(|(i, is): (usize, &Vec<TLweSample>)| {
        is.iter()
          .map(|js: &TLweSample| {
            let new_a: Vec<TorusPolynomial> = js
              .a
              .iter()
              .map(|a: &TorusPolynomial| {
                let new_coefs = a
                  .coefs
                  .iter()
                  .enumerate()
                  .map(
                    |(coef_idx, coef): (usize, &i32)| {
                      if coef_idx == 0 {
                        coef + h[i]
                      } else {
                        *coef
                      }
                    },
                  )
                  .collect::<Vec<Torus32>>();
                TorusPolynomial { coefs: new_coefs }
              })
              .collect();
            TLweSample {
              a: new_a,
              ..js.clone()
            }
          })
          .collect()
      })
      .collect();

    // Is equivalent to:

    // for i in 0..l as usize {
    //   for j in 0..=k as usize {
    //     self.all_sample[i][j].a[j].coefs[0] += message * h[i];
    //   }
    // }
  }

  #[allow(clippy::needless_range_loop)]
  pub(crate) fn add_h(&mut self, params: &TGswParams) {
    let k = params.tlwe_params.k;
    let l = params.l;
    let h = &params.h;

    // compute self += H
    for i in 0..l as usize {
      for j in 0..=k as usize {
        self.all_sample[i][j].a[j].coefs[0] += h[i];
      }
    }
  }

  #[allow(clippy::needless_range_loop)]
  pub(crate) fn add_mu_h(&mut self, message: &IntPolynomial, params: &TGswParams) {
    let k = params.tlwe_params.k;
    let n = params.tlwe_params.n;
    let l = params.l;
    let h = &params.h;
    let mu = &message.coefs;

    // Compute self += H
    for i in 0..l as usize {
      for j in 0..=k as usize {
        let target = &mut self.all_sample[i][j].a[j].coefs;
        println!("target coefs befor: {:?}", target);
        target
          .iter_mut()
          .zip(mu.iter())
          .for_each(|(t, mu)| *t += mu * h[i]);
        println!("target coefs after: {:?}", target);

        // for jj in 0..n as usize {
        //   println!(
        //     "Target len: {}, mu len: {}, h len: {}, jj: {}, n: {}",
        //     target.len(),
        //     mu.len(),
        //     h.len(),
        //     jj,
        //     n
        //   );
        //   target[jj] += mu[jj] * h[i];
        // }
      }
    }
  }
}

/// Update l'accumulateur ligne 5 de l'algo toujours
/// void tGswTLweDecompH(IntPolynomial* result, const TLweSample* sample,const TGswParams* params);
/// accum *= sample
pub(crate) fn tgsw_extern_mul_to_tlwe(
  accum: &mut TLweSample,
  sample: &TGswSample,
  params: &TGswParams,
) {
  let par = &params.tlwe_params;
  let n = par.n;
  let kpl = params.kpl;

  // TODO: improve this new/delete
  //   IntPolynomial *dec = new_IntPolynomial_array(kpl, N);

  let mut dec = vec![vec![IntPolynomial::new(n); params.l as usize]; (par.k + 1) as usize];
  tgsw_tlwe_decomposition_h(&mut dec, accum, params);

  // TODO: Remove this and remove mutability-requiring functions
  accum.clear();
  dec
    .iter()
    .zip(sample.all_sample.iter())
    .for_each(|(d, a)| accum.add_mul_r_(&d[0], &a[0], par));
  // for i in 0..dec.len() as usize {
  //   println!("kpl: {}, k: {}, l: {}, i: {}", kpl, par.k, params.l, i);
  //   // TODO: Figure out if this is supposed to be [0][i] instead, or something else...
  //   let d = &dec[i][0];
  //   let ass = &sample.all_sample[i][0];
  //   accum.add_mul_r_(d, ass, par);
  // }

  //   for (int32_t i = 0; i < kpl; i++) {
  //       tLweAddMulRTo(accum, &dec[i], &sample->all_sample[i], par);
  //   }
}

/// Fonction de decomposition
fn tgsw_tlwe_decomposition_h(
  result: &mut Vec<Vec<IntPolynomial>>,
  sample: &mut TLweSample,
  params: &TGswParams,
) {
  let k = params.tlwe_params.k;
  let l = params.l;
  for i in 0..=k {
    // b=a[k]
    tgsw_torus32_polynomial_decomposition_h(
      &mut result[(i/* /* TODO: Remove this when you figure this out: Don't think this is necessary? */ * l*/)
        as usize],
      &mut sample.a[i as usize],
      params,
    );
    //     tGswTorus32PolynomialDecompH(result + (i * l), &sample->a[i], params);
  }
}

fn tgsw_torus32_polynomial_decomposition_h(
  result: &mut Vec<IntPolynomial>,
  sample: &mut TorusPolynomial,
  params: &TGswParams,
) {
  let n = params.tlwe_params.n;
  let l = params.l;
  let bg_bit = params.bg_bit;
  let buf = &mut sample.coefs;
  let mask_mod = params.mask_mod;
  let half_bg = params.half_bg;
  let offset = params.offset;
  // First, add offset to everyone
  for j in 0..n as usize {
    buf[j] += offset as i32;
  }

  // Then, do the decomposition (TODO: in parallel)
  for p in 0..l as usize {
    let decal = 32 - (p + 1) as i32 * bg_bit;
    let res_p = &mut result[p].coefs;
    for j in 0..n as usize {
      let temp1 = (buf[j] >> decal) & mask_mod as i32;
      res_p[j] = temp1 - half_bg;
    }
  }
}

pub struct TGswSampleFFT {
  /// TLweSample* all_sample; (k+1)l TLwe Sample
  all_samples: Vec<TLweSampleFFT>,
  //   TLweSampleFFT **sample; /// optional access to the different size blocks l
  sample: Option<TLweSampleFFT>,
  k: i32,
  l: i32,
}

// impl TGswSampleFFT {
//   fn new(params: TGswParams, all_samples_raw: Vec<TLweSampleFFT>) -> Self {
//     let k = params.tlwe_params.k;
//     let l = params.l;
//     let all_samples = all_samples_raw;
//     let mut sample: Vec<Option<TLweSampleFFT>> = vec![None; ((k + 1) * l) as usize];
//     for p in 0..(k + 1) {
//       sample[p as usize] = Some(all_samples[(p * l) as usize]);
//     }
//     Self {
//       all_samples,
//       sample,
//       k,
//       l,
//     }
//   }
// }

#[cfg(test)]
mod tests {
  use super::*;
  use crate::tlwe::TorusPolynomial;

  fn generate_parameters() -> Vec<TGswParams> {
    vec![
      TGswParams::new(4, 8, TLweParameters::new(512, 1, 0f64, 1f64)),
      TGswParams::new(3, 10, TLweParameters::new(512, 2, 0f64, 1f64)),
      TGswParams::new(3, 10, TLweParameters::new(1024, 1, 0f64, 1f64)),
      TGswParams::new(4, 8, TLweParameters::new(1024, 2, 0f64, 1f64)),
      TGswParams::new(4, 8, TLweParameters::new(2048, 1, 0f64, 1f64)),
      TGswParams::new(3, 10, TLweParameters::new(2048, 2, 0f64, 1f64)),
    ]
  }

  fn generate_keys() -> Vec<TGswKey> {
    generate_parameters().iter().map(TGswKey::new).collect()
  }

  /*
   *  Testing the function tGswKeyGen
   * This function generates a random TLwe key for the given parameters
   * The TLwe key for the result must be allocated and initialized
   * (this means that the parameters are already in the result)
   */
  #[test]
  fn test_key_generation() {
    for param in generate_parameters() {
      let mut key = TGswKey::new(&param);
      key.tlwe_key.generate();

      // Assert key is binary (could be eliminated by )
      assert!(key
        .key
        .iter()
        .all(|key| key.coefs.iter().all(|&c| c == 1i32 || c == 0i32)));
    }
  }

  fn fully_random_tgsw(sample: &mut TGswSample, alpha: f64, params: &TGswParams) {
    let l = params.l;
    let k = params.tlwe_params.k;

    // This is butt-ugly
    for j in 0..l as usize {
      for i in 0..=k as usize {
        let mut row = &mut sample.all_sample[j][i];
        for u in 0..=k {
          row.a[u as usize] = TorusPolynomial::torus_polynomial_uniform(row.a.len() as i32);
        }
        row.current_variance = alpha * alpha;
      }
    }
  }

  #[test]
  fn test_add_h() {
    for params in generate_parameters() {
      let mut sample = TGswSample::new(&params);
      let kpl = params.kpl;
      let l = params.l;
      let k = params.tlwe_params.k;
      let n = params.tlwe_params.n;
      let h = &params.h;
      let alpha = 4.2;

      fully_random_tgsw(&mut sample, alpha, &params);

      let sample_copy = sample.clone();
      sample.add_h(&params);

      // Verify all coefficients
      for i in 0..l as usize {
        for j in 0..=k as usize {
          assert_eq!(
            sample.all_sample[i][j].current_variance,
            sample_copy.all_sample[i][j].current_variance
          );

          for u in 0..=k as usize {
            //verify that pol[bloc][i][u]=initial[bloc][i][u]+(bloc==u?hi:0)
            let new_polynomial = &sample.all_sample[i][j].a[u];
            let old_polynomial = &sample_copy.all_sample[i][j].a[u];
            assert_eq!(
              new_polynomial.coefs[0], // Should this be i == u?
              old_polynomial.coefs[0] + (if j == u { h[i] } else { 0 })
            );
            assert_eq!(new_polynomial.coefs[1..], old_polynomial.coefs[1..]);
          }
        }
      }
    }
  }

  fn random_int_polynomial(n: i32) -> IntPolynomial {
    let mut rng = rand::thread_rng();
    use rand::distributions::Distribution;
    let d = rand_distr::Uniform::new(i32::min_value(), i32::max_value());

    let coefs: Vec<i32> = (0..n).map(|_| d.sample(&mut rng) % 10 - 5).collect();
    assert_eq!(coefs.len() as i32, n);
    IntPolynomial { coefs, n }
  }

  #[test]
  #[ignore]
  fn test_add_mu_h() {
    for params in generate_parameters() {
      let mut sample = TGswSample::new(&params);
      let kpl = params.kpl;
      let l = params.l;
      let k = params.tlwe_params.k;
      let n = params.tlwe_params.n;
      let h = &params.h;
      let alpha = 4.2;
      let message = random_int_polynomial(n);

      fully_random_tgsw(&mut sample, alpha, &params);

      let sample_copy = sample.clone();
      sample.add_mu_h(&message, &params);

      // Verify all coefficients
      for i in 0..l as usize {
        for j in 0..=k as usize {
          assert_eq!(
            sample.all_sample[i][j].current_variance,
            sample_copy.all_sample[i][j].current_variance
          );

          for u in 0..=k as usize {
            //verify that pol[bloc][i][u]=initial[bloc][i][u]+(bloc==u?hi*mess:0)
            let new_polynomial = &sample.all_sample[i][j].a[u];
            let old_polynomial = &sample_copy.all_sample[i][j].a[u];

            if j == u {
              new_polynomial
                .coefs
                .iter()
                .zip(old_polynomial.coefs.iter())
                .zip(message.coefs.iter())
                .for_each(|((n, o), m)| assert_eq!(*n, *o + h[i] * (dbg!(*m))));
            // for jj in 0..n as usize {
            //   assert_eq!(
            //     new_polynomial.coefs[jj],
            //     old_polynomial.coefs[jj] + h[i] * message.coefs[jj]
            //   );
            // }
            } else {
              assert!(new_polynomial
                .coefs
                .iter()
                .zip(old_polynomial.coefs.iter())
                .all(|(a, b)| a == b));
            }
            assert_eq!(
              new_polynomial.coefs[0], // Should this be i == u?
              old_polynomial.coefs[0] + (if j == u { dbg!(h[i]) } else { 0 })
            );
            assert_eq!(new_polynomial.coefs[1..], old_polynomial.coefs[1..]);
          }
        }
      }
    }
  }
}
