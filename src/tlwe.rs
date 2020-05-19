use crate::lwe::{LweParams, LweSample};
use crate::numerics::{gaussian32, torus_polynomial_mul_r};
use crate::polynomial::{IntPolynomial, Polynomial, TorusPolynomial};

use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
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
    use rand::Rng;
    let mut rng = rand::thread_rng();

    // Fill key with random integers
    let key: Vec<IntPolynomial> = (0..k)
      .map(|_| {
        IntPolynomial::from(
          (0..n)
            .map(|_| if rng.gen::<bool>() { 1 } else { 0 })
            .collect::<Vec<i32>>(),
        )
      })
      .collect();

    Self {
      key,
      params: params.clone(),
    }
  }
}

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
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
      current_variance: 0_f64,
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
    let a_last_part =
      TorusPolynomial::from((0..n).map(|_| gaussian32(0, alpha)).collect::<Vec<i32>>());

    // Multiply key.key with self.a and sum them up, add it to b (last part)
    let poly_sum = key
      .key
      .iter()
      .zip(self.a.iter())
      .map(|(a, b)| crate::numerics::poly_multiplier(a, b))
      .fold(TorusPolynomial::zero(n as usize), |acc, p| acc + p);

    self.a = a_part
      .chain(std::iter::once(a_last_part + poly_sum))
      .collect();

    self.current_variance = alpha * alpha;
  }

  /// Sets all values to zero
  /// TODO: Remove this function when things stabilize. We shouldn't need to reuse values like in C++/C
  pub(crate) fn clear(&mut self) {
    self.a = self
      .a
      .iter()
      .map(|poly| TorusPolynomial::zero(poly.len()))
      .collect();

    self.current_variance = 0_f64;
  }

  /// self += p * sample
  pub(crate) fn add_mul_r_(&mut self, p: &IntPolynomial, sample: &Self, params: &TLweParameters) {
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
      l.coefficients[(i * n) as usize] = self.a[i as usize][0];
      for j in 1..n {
        // l->a[i*N+j] = -x->a[i].coefsT[N+0-j];
        l.coefficients[(i * n + j) as usize] = -self.a[i as usize][(n - j) as usize];
      }
    }

    match self.a.len() {
      0 => panic!("Cannot get last element of a, as it is empty!"),
      n => l.b = self.a[n - 1][0],
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

pub(crate) fn mul_by_monomial(x: TLweSample, shift: i32) -> TLweSample {
  TLweSample {
    a: x
      .a
      .into_iter()
      .map(|c| crate::polynomial::mul_by_monomial(c, shift))
      .collect(),
    current_variance: x.current_variance,
    k: x.k,
  }
}
