use crate::tlwe::{IntPolynomial, TLweKey, TLweParameters, TLweSample, TLweSampleFFT, Torus32};

#[derive(Clone)]
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
  tlwe_params: TLweParameters,

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
  params: TGswParams,

  /// the tlwe params of each rows
  tlwe_params: TLweParameters,

  /// the key (array of k polynomials)
  key: Vec<IntPolynomial>,

  tlwe_key: TLweKey,
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

  pub(crate) fn encrypt(&mut self, result: &mut TGswSample, message: i32, alpha: f64) {
    result.encrypt_zero(alpha, &self);
    result.add_mu_int_h(message, &self.params);
  }
}

#[derive(Clone, Debug)]
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
    let l = params.l; // Lines / rows
                      // TODO: find out if this is correctamente
    let all_sample = vec![vec![TLweSample::new(&params.tlwe_params); (k + 1) as usize]; l as usize];

    Self { all_sample, k, l }
  }

  pub(crate) fn encrypt_zero(&mut self, alpha: f64, key: &TGswKey) {
    let rl_key = &key.tlwe_key;
    let kpl = key.params.kpl;
    for p in 0..kpl {
      // self.all_sample[p as usize].encrypt_zero(alpha, rl_key);
    }
    unimplemented!()
  }

  pub(crate) fn add_mu_int_h(&mut self, message: i32, params: &TGswParams) {
    let k = params.tlwe_params.k;
    let l = params.l;
    let h = &params.h;

    // Compute self += H
    for i in 0..l as usize {
      for j in 0..=k as usize {
        self.all_sample[i][j].a[j].coefs[0] += message * h[i];
      }
    }
  }

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
}
