use crate::tlwe::{IntPolynomial, TLweKey, TLweParameters, TLweSample, TLweSampleFFT, Torus32};

struct TGswParams {
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
  fn new(l: i32, bg_bit: i32, tlwe_params: TLweParameters) -> Self {
    let bg = (1 << bg_bit) as i32;
    let half_bg = (bg >> 1) as i32;
    let mask_mod = (bg - 1) as u32;
    let mut h = vec![l];
    for i in 0..l {
      let kk = (32 - (i + 1)) * bg_bit;
      // 1/(Bg^(i+1)) as a Torus32
      h.push(1 << kk);
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

struct TGswKey {
  /// the parameters
  params: TGswParams,

  /// the tlwe params of each rows
  tlwe_params: TLweParameters,

  /// the key (array of k polynomials)
  key: IntPolynomial,

  tlwe_key: TLweKey,
}

// impl TGswKey {
//   // same key as in TLwe
//   fn new(params: TGswParams) -> Self {
//     let tlwe_params = params.tlwe_params;
//     let tlwe_key = tlwe_params;
//     let key = tlwe_key.key;
//     Self {
//       params,
//       tlwe_params,
//       key,
//       tlwe_key,
//     }
//   }
// }

struct TGswSample {
  /// (k+1)l TLwe Sample
  all_sample: TLweSample,
  /// optional access to the different size blocks l
  bloc_sample: Vec<Vec<TLweSample>>,
  current_variance: f32,
  k: i32,
  l: i32,
}

struct TGswSampleFFT {
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
