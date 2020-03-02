use crate::numerics::gaussian32;
use rand::distributions::Distribution;
#[derive(Clone)]
pub struct LWEParams {
  n: i32,
  alpha_min: f64,
  alpha_max: f64,
}

impl LWEParams {
  pub fn new(n: i32, alpha_min: f64, alpha_max: f64) -> Self {
    Self {
      n,
      alpha_min,
      alpha_max,
    }
  }
}

#[derive(Clone)]
pub struct TLweParameters {
  /// a power of 2: degree of the polynomials
  pub n: i32,
  /// number of polynomials in the mask
  pub k: i32,
  /// minimal noise s.t. the sample is secure
  pub alpha_min: f64,
  /// maximal noise s.t. we can decrypt
  pub alpha_max: f64,
  /// lwe params if one extracts
  pub extracted_lweparams: LWEParams,
}

impl TLweParameters {
  pub fn new(n: i32, k: i32, alpha_min: f64, alpha_max: f64) -> Self {
    Self {
      n,
      k,
      alpha_min,
      alpha_max,
      extracted_lweparams: LWEParams::new(n * k, alpha_min, alpha_max),
    }
  }
}

#[derive(Clone)]
pub struct IntPolynomial {
  pub(crate) n: i32,
  pub(crate) coefs: Vec<i32>,
}

impl IntPolynomial {
  pub(crate) fn new(n: i32) -> Self {
    Self {
      n,
      coefs: vec![0; n as usize],
    }
  }
}

pub struct TLweKey {
  /// Parameters of the key
  params: TLweParameters,
  /// the key (i.e k binary polynomials)
  pub(crate) key: Vec<IntPolynomial>,
}

impl TLweKey {
  pub(crate) fn new(params: &TLweParameters) -> Self {
    Self {
      params: params.clone(),
      key: vec![IntPolynomial::new(params.n); params.n as usize],
    }
  }

  pub(crate) fn generate(&mut self) {
    let n = self.params.n;
    let k = self.params.k;
    let d = rand_distr::Uniform::new(0, 1);
    let mut rng = rand::thread_rng();
    for i in 0..k {
      for j in 0..n {
        self.key[i as usize].coefs[j as usize] = d.sample(&mut rng);
      }
    }
  }
}

/// Idea:
/// we may want to represent an element x of the real torus by
/// the integer rint(2^32.x) modulo 2^32
///  -- addition, subtraction and integer combinations are native operation
///  -- modulo 1 is mapped to mod 2^32, which is also native!
/// This looks much better than using float/doubles, where modulo 1 is not
/// natural at all.
pub type Torus32 = i32;

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct TorusPolynomial {
  pub(crate) n: i32,
  pub(crate) coefs: Vec<Torus32>,
}

impl TorusPolynomial {
  pub(crate) fn new(n: i32) -> Self {
    Self {
      n,
      coefs: vec![0; n as usize],
    }
  }

  pub(crate) fn torus_polynomial_uniform(n: i32) -> Self {
    let d = rand_distr::Uniform::new(i32::min_value(), i32::max_value());
    let mut rng = rand::thread_rng();

    Self {
      n,
      coefs: (0..n as i32).map(|_| d.sample(&mut rng)).collect(),
    }
  }
}

#[derive(Clone, Debug, PartialEq)]
pub struct TLweSample {
  /// array of length k+1: mask + right term
  pub(crate) a: Vec<TorusPolynomial>,
  /// FIXME: This is some C++ shit. b is actually referring to a single value within a
  /// alias of a[k] to get the right term
  b: TorusPolynomial,
  /// avg variance of the sample
  pub(crate) current_variance: f64,
  k: i32,
}
impl TLweSample {
  pub(crate) fn new(params: &TLweParameters) -> Self {
    //Small change here:
    //a is a table of k+1 polynomials, b is an alias for &a[k]
    //like that, we can access all the coefficients as before:
    //  &sample->a[0],...,&sample->a[k-1]  and &sample->b
    //or we can also do it in a single for loop
    //  &sample->a[0],...,&sample->a[k]
    let a = vec![TorusPolynomial::new(params.n); (params.k + 1) as usize];
    let b = a[params.k as usize].clone();
    Self {
      a,
      b,
      current_variance: 0f64,
      k: params.k,
    }
  }

  /// Create an homogeneous tlwe sample
  pub(crate) fn encrypt_zero(&mut self, alpha: f64, key: &TLweKey) {
    let n = key.params.n;
    let k = key.params.k;

    self.b.coefs = vec![gaussian32(0, alpha); n as usize];

    // Random-generate tori
    self.a = (0..k)
      .map(|_| TorusPolynomial::torus_polynomial_uniform(n))
      .collect();

    // torusPolynomialAddMulR(result->b, &key->key[i], &result->a[i]);

    for i in 0..k {
      crate::numerics::torus_polynomial_mul_r(
        &mut self.b,
        &key.key[i as usize],
        &self.a[i as usize],
      );
    }

    self.current_variance = alpha * alpha;

    unimplemented!()
  }
}

/// Figure out what this is
#[derive(Clone)]
struct LagrangeHalfCPolynomial<Data, Precomp> {
  data: Data,
  precomp: Precomp,
}

#[derive(Clone)]
pub struct TLweSampleFFT {
  /// array of length k+1: mask + right term
  a: Vec<LagrangeHalfCPolynomial<u8, u8>>,
  /// FIXME: This is some C++ shit. b is actually referring to a single value within a
  /// alias of a[k] to get the right term
  b: LagrangeHalfCPolynomial<u8, u8>,
  /// avg variance of the sample
  current_variance: f64,
  /// TODO: Figure out if this is required...
  /// required during the destructor call...
  k: i32,
}
impl TLweSampleFFT {
  fn new(
    params: TLweParameters,
    arr: Vec<LagrangeHalfCPolynomial<u8, u8>>,
    current_variance: f64,
  ) -> Self {
    let k = params.k;
    let b = (&arr[k as usize]).clone();
    Self {
      a: arr,
      // a is a table of k+1 polynomials, b is an alias for &a[k]
      b,
      current_variance,
      k,
    }
  }
}
