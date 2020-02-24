#[derive(Clone)]
pub(crate) struct LWEParams {
  n: i32,
  alpha_min: f32,
  alpha_max: f32,
}

#[derive(Clone)]
pub(crate) struct TLweParameters {
  /// a power of 2: degree of the polynomials
  pub(crate) n: i32,
  /// number of polynomials in the mask
  pub(crate) k: i32,
  /// minimal noise s.t. the sample is secure
  pub(crate) alpha_min: f32,
  /// maximal noise s.t. we can decrypt
  pub(crate) alpha_max: f32,
  /// lwe params if one extracts
  pub(crate) extracted_lweparams: LWEParams,
}

impl TLweParameters {
  fn new(n: i32, k: i32, alpha_min: f32, alpha_max: f32) -> Self {
    Self {
      n,
      k,
      alpha_min,
      alpha_max,
      extracted_lweparams: LWEParams {
        n: n * k,
        alpha_min,
        alpha_max,
      },
    }
  }
}

#[derive(Clone)]
pub(crate) struct IntPolynomial {
  n: i32,
  coefs: Vec<i32>,
}

impl IntPolynomial {
  fn new(n: i32) -> Self {
    Self {
      n,
      coefs: vec![0; n as usize],
    }
  }
}

pub(crate) struct TLweKey {
  /// Parameters of the key
  params: TLweParameters,
  /// the key (i.e k binary polynomials)
  key: Vec<IntPolynomial>,
}

impl TLweKey {
  fn new(params: TLweParameters) -> Self {
    Self {
      params: params.clone(),
      key: vec![IntPolynomial::new(params.n); params.n as usize],
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
pub(crate) type Torus32 = i32;

struct TorusPolynomial {
  n: i32,
  coefs: Vec<Torus32>,
}

pub(crate) struct TLweSample {
  /// array of length k+1: mask + right term
  a: TorusPolynomial,
  /// FIXME: This is some C++ shit. b is actually referring to a single value within a
  /// alias of a[k] to get the right term
  b: TorusPolynomial,
  /// avg variance of the sample
  current_variance: f32,
  k: i32,
}

// impl TLweSample {
//   fn new(params: TLweParameters) -> Self {
//     //Small change here:
//     //a is a table of k+1 polynomials, b is an alias for &a[k]
//     //like that, we can access all the coefficients as before:
//     //  &sample->a[0],...,&sample->a[k-1]  and &sample->b
//     //or we can also do it in a single for loop
//     //  &sample->a[0],...,&sample->a[k]
//     Self {

//     }
//   }
// }

/// Figure out what this is
#[derive(Clone)]
struct LagrangeHalfCPolynomial<Data, Precomp> {
  data: Data,
  precomp: Precomp,
}

#[derive(Clone)]
pub(crate) struct TLweSampleFFT {
  /// array of length k+1: mask + right term
  a: Vec<LagrangeHalfCPolynomial<u8, u8>>,
  /// FIXME: This is some C++ shit. b is actually referring to a single value within a
  /// alias of a[k] to get the right term
  b: LagrangeHalfCPolynomial<u8, u8>,
  /// avg variance of the sample
  current_variance: f32,
  /// TODO: Figure out if this is required...
  /// required during the destructor call...
  k: i32,
}
impl TLweSampleFFT {
  fn new(
    params: TLweParameters,
    arr: Vec<LagrangeHalfCPolynomial<u8, u8>>,
    current_variance: f32,
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
