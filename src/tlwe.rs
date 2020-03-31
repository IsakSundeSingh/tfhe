use crate::lwe::{LweParams, LweSample};
use crate::numerics::{
  gaussian32, torus_polynomial_mul_by_xai_minus_one, torus_polynomial_mul_r, Modulo, Torus32,
};
use crate::polynomial::{IntPolynomial, Polynomial, TorusPolynomial};
use rand::distributions::Distribution;

#[derive(Clone, Debug, PartialEq)]
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
  pub extracted_lweparams: LweParams,
}

impl TLweParameters {
  pub fn new(n: i32, k: i32, alpha_min: f64, alpha_max: f64) -> Self {
    Self {
      n,
      k,
      alpha_min,
      alpha_max,
      extracted_lweparams: LweParams::new(n * k, alpha_min, alpha_max),
    }
  }
}

pub struct TLweKey {
  /// Parameters of the key
  pub(crate) params: TLweParameters,
  /// the key (i.e k binary polynomials)
  pub(crate) key: Vec<IntPolynomial>,
}

impl TLweKey {
  /// Generates a random key with the specified parameters.
  pub(crate) fn generate(params: &TLweParameters) -> Self {
    let n = params.n;
    let k = params.k;
    let d = rand_distr::Uniform::new(0, 1);
    let mut rng = rand::thread_rng();

    // Fill key with random integers
    let key: Vec<IntPolynomial> = (0..k)
      .map(|_| IntPolynomial::from((0..n).map(|_| d.sample(&mut rng)).collect::<Vec<i32>>()))
      .collect();

    Self {
      key,
      params: params.clone(),
    }
  }
}

#[derive(Clone, Debug, PartialEq)]
pub struct TLweSample {
  /// array of length k+1: mask + right term
  pub(crate) a: Vec<TorusPolynomial>,
  // This field was actually the last element of `a`
  // pub(crate) b: TorusPolynomial,
  /// avg variance of the sample
  pub(crate) current_variance: f64,
  pub(crate) k: i32,
}

impl TLweSample {
  pub(crate) fn new(params: &TLweParameters) -> Self {
    //Small change here:
    //a is a table of k+1 polynomials, b is an alias for &a[k]
    //like that, we can access all the coefficients as before:
    //  &sample->a[0],...,&sample->a[k-1]  and &sample->b
    //or we can also do it in a single for loop
    //  &sample->a[0],...,&sample->a[k]
    let a = vec![TorusPolynomial::new(params.n); (params.k + 2) as usize];
    Self {
      a,
      current_variance: 0f64,
      k: params.k,
    }
  }

  /// Creates a noiseless `TLweSample` with a given μ
  ///
  /// # Panics
  ///
  /// Panics if the number of elements in the new sample is zero
  pub(crate) fn trivial(mu: TorusPolynomial, params: &TLweParameters) -> Self {
    let mut sample = Self::new(params);
    match sample.a.len() {
      // If the length is empty, panic
      0 => {
        panic!("Trying to create a new TLweSample with a mu value, but number of elements is zero!")
      }
      n => sample.a[n - 1] = mu,
    }
    sample
  }

  /// Create an homogeneous tlwe sample
  pub(crate) fn encrypt_zero(&mut self, alpha: f64, key: &TLweKey) {
    let n = key.params.n;
    let k = key.params.k;

    // Random-generate tori
    let a_part = (0..k).map(|_| TorusPolynomial::uniform(n as usize));
    let a_last_part = TorusPolynomial::from(vec![gaussian32(0, alpha); n as usize]);

    // Multiply key.key with self.a and sum them up, add it to b (last part)
    let poly_sum = key
      .key
      .iter()
      .zip(self.a.iter())
      .map(|(a, b)| crate::numerics::poly_multiplier(a, b))
      .fold(TorusPolynomial::zero(n as usize), |acc, p| acc + p);

    let a = a_part
      .chain(std::iter::once(a_last_part + poly_sum))
      .collect();
    self.a = a;

    // for i in 0..k {
    //   crate::numerics::torus_polynomial_mul_r(
    //     &mut self.b,
    //     &key.key[i as usize],
    //     &self.a[i as usize],
    //   );
    // }

    self.current_variance = alpha * alpha;
  }

  fn encrypt_sym_t(&mut self, message: Torus32, alpha: f64, key: &TLweKey) {
    self.encrypt_zero(alpha, key);
    match self.a.len() {
      0 => panic!("Could not set b as a had length 0!"),
      n => self.a[n - 1].coefs[0] += message,
    }
  }

  /// Sets all values to zero
  /// TODO: Remove this function when things stabilize. We shouldn't need to reuse values like in C++/C
  pub(crate) fn clear(&mut self) {
    self.a = self
      .a
      .iter()
      .map(|poly| TorusPolynomial::zero(poly.len()))
      .collect();

    self.current_variance = 0f64;
  }

  /// self += p * sample
  pub(crate) fn add_mul_r_(
    &mut self,
    p: &IntPolynomial,
    sample: &TLweSample,
    params: &TLweParameters,
  ) {
    let k = params.k;

    for i in 0..=k as usize {
      torus_polynomial_mul_r(&mut self.a[i], p, &sample.a[i]);
    }

    self.current_variance += p.norm_squared() * sample.current_variance;
  }

  pub(crate) fn extract_lwe(self, params: &LweParams, rparams: &TLweParameters) -> LweSample {
    let n = rparams.n;
    let k = rparams.k;

    // TODO: This might be incorrect and unnecessary
    assert_eq!(params.n, k * n);
    let mut l = LweSample::new(params);

    for i in 0..k {
      l.coefficients[(i * n) as usize] = self.a[i as usize].coefs()[0];
      for j in 1..n {
        // l->a[i*N+j] = -x->a[i].coefsT[N+0-j];
        l.coefficients[(i * n + j) as usize] = -self.a[i as usize].coefs()[(n - j) as usize];
      }
    }

    match self.a.len() {
      0 => panic!("Cannot get last element of a, as it is empty!"),
      n => l.b = self.a[n - 1].coefs()[0],
    }
    l
  }
}

impl std::ops::Add<TLweSample> for TLweSample {
  type Output = Self;
  fn add(self, sample: Self) -> Self {
    Self {
      a: self
        .a
        .into_iter()
        .zip(sample.a.into_iter())
        .map(|(a, b)| a + b)
        .collect(),
      current_variance: self.current_variance + sample.current_variance,
      ..self
    }
  }
}

impl std::ops::Sub<TLweSample> for TLweSample {
  type Output = Self;
  #[allow(clippy::suspicious_arithmetic_impl)]
  fn sub(self, sample: Self) -> Self {
    Self {
      a: self
        .a
        .into_iter()
        .zip(sample.a.into_iter())
        .map(|(a, b)| a - b)
        .collect(),
      current_variance: self.current_variance + sample.current_variance,
      ..self
    }
  }
}

/// TODO: Remove this as it is ugly
pub(crate) fn tlwe_add_to(res: &mut TLweSample, sample: &TLweSample) {
  *res = TLweSample {
    a: res
      .a
      .iter()
      .zip(sample.a.iter())
      // TODO: Fix inefficient memory allocations (clone)
      .map(|(a, b): (&TorusPolynomial, &TorusPolynomial)| a.clone() + b.clone())
      .collect(),
    current_variance: res.current_variance + sample.current_variance,
    k: res.k,
  }
}

pub(crate) fn mul_by_monomial(x: TLweSample, shift: i32) -> TLweSample {
  TLweSample {
    a: x
      .a
      .into_iter()
      .map(|c| {
        TorusPolynomial::from(crate::polynomial::mul_by_monomial(
          IntPolynomial::from(c),
          shift,
        ))
      })
      .collect(),
    current_variance: x.current_variance,
    k: x.k,
  }
}

/// Mult externe de X^ai-1 par bki
pub(crate) fn tlwe_mul_by_xai_minus_one(
  result: &TLweSample,
  ai: i32,
  bk: &TLweSample,
  params: &TLweParameters,
) -> TLweSample {
  let torus_polynomials = bk
    .a
    .iter()
    .map(|ba| torus_polynomial_mul_by_xai_minus_one(ai.modulo(2 * params.n), ba))
    .collect();

  // for i in 0..=k as usize {
  //   torus_polynomial_mul_by_xai_minus_one(ai, &bk.a[i]);
  //   // torusPolynomialMulByXaiMinusOne(&result->a[i], ai, &bk->a[i]);
  // }

  TLweSample {
    a: torus_polynomials,
    ..result.clone()
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
