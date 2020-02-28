use crate::tlwe::TorusPolynomial;
use crate::tlwe::{IntPolynomial, Torus32};

// Gaussian sample centered in message, with standard deviation sigma
pub(crate) fn gaussian32(message: Torus32, sigma: f64) -> Torus32 {
  use rand::distributions::Distribution;

  // Attention: all the implementation will use the stdev instead of the gaussian fourier param
  let d = rand_distr::Normal::new(0f64, sigma).expect("Could not create normal distribution");
  let mut rng = rand::thread_rng();
  let error = d.sample(&mut rng);

  // Overflowed here, using wrapping add to imitate C++ behavior
  message.wrapping_add(f64_to_torus_32(error))
}

/// # Warning
/// Weird and lossy conversion
pub(crate) fn f64_to_torus_32(d: f64) -> Torus32 {
  // 2 ^ 32
  const TWO_32: i64 = 4_294_967_296;
  let inner = d - ((d as i64) as f64);
  let x = (inner * TWO_32 as f64) as i64;
  x as Torus32
  // C++ conversion: // return int32_t(int64_t((d - int64_t(d))*_two32));
}

pub(crate) fn torus_32_to_f64(x: Torus32) -> f64 {
  // 2 ^ 32
  const TWO_32: f64 = 4_294_967_296.0;
  (x as f64) / TWO_32
}

/// Used to approximate the phase to the nearest message possible in the message space
/// The constant Msize will indicate on which message space we are working (how many messages possible)
///
/// "travailler sur 63 bits au lieu de 64, car dans nos cas pratiques, c'est plus précis"
pub(crate) fn approximate_phase(phase: Torus32, message_size: i32) -> Torus32 {
  // Width of each interval
  let interval: u64 = ((1u64 << 63) / (message_size as u64)) * 2;

  // Beginning of the first interval
  let half_interval: u64 = interval / 2;

  // Overflowed here, using wrapping add to imitate C++ behavior
  let mut phase64: u64 = ((phase as u64) << 32).wrapping_add(half_interval);

  // Floor to the nearest multiples of interval
  phase64 -= phase64 % interval;

  // Rescale to torus32
  (phase64 >> 32) as Torus32
}

/// Used to approximate the phase to the nearest message possible in the message space
/// The constant Msize will indicate on which message space we are working (how many messages possible)
///
/// "travailler sur 63 bits au lieu de 64, car dans nos cas pratiques, c'est plus précis"
pub(crate) fn mod_switch_to_torus32(mu: i32, message_size: i32) -> Torus32 {
  // Width of each interval
  let interval: u64 = ((1u64 << 63) / message_size as u64) * 2;

  // Overflowed here, using wrapping mul to imitate C++ behavior
  let phase64: u64 = (mu as u64).wrapping_mul(interval);

  // Floor to the nearest multiples of interval
  (phase64 >> 32) as Torus32
}

pub(crate) fn torus_polynomial_mul_r(
  result: &mut TorusPolynomial,
  poly1: &IntPolynomial,
  poly2: &TorusPolynomial,
) {
  let res = poly_multiplier(poly1, poly2);
  assert_eq!(result.n, res.n);

  result.coefs = result
    .coefs
    .iter()
    .zip(res.coefs.iter())
    .map(|(a, b)| a + b)
    .collect();
  // let tmp = crate::tlwe::LagrangeHalfCPolynomial
  // const int32_t N = poly1->N;
  // LagrangeHalfCPolynomial* tmp = new_LagrangeHalfCPolynomial_array(3,N);
  // TorusPolynomial* tmpr = new_TorusPolynomial(N);
  // IntPolynomial_ifft(tmp+0,poly1);
  // TorusPolynomial_ifft(tmp+1,poly2);
  // LagrangeHalfCPolynomialMul(tmp+2,tmp+0,tmp+1);
  // TorusPolynomial_fft(tmpr, tmp+2);
  // torusPolynomialAddTo(result, tmpr);
  // delete_TorusPolynomial(tmpr);
  // delete_LagrangeHalfCPolynomial_array(3,tmp);
}

fn poly_multiplier(a: &IntPolynomial, b: &TorusPolynomial) -> TorusPolynomial {
  assert_eq!(a.n, a.coefs.len() as i32);
  assert_eq!(b.n, b.coefs.len() as i32);

  let degree = a.n + b.n - 2;
  let mut coefs = vec![0; (degree + 1) as usize];

  for i in 0..a.n {
    for j in 0..b.n {
      coefs[(i + j) as usize] += a.coefs[i as usize] * b.coefs[j as usize];
    }
  }

  TorusPolynomial {
    n: coefs.len() as i32,
    coefs,
  }
}

#[test]
fn test_poly_multiplier() {
  let a = IntPolynomial {
    n: 3,
    coefs: vec![10, 20, 30],
  };

  let b = TorusPolynomial {
    n: 3,
    coefs: vec![1, 2, 3],
  };

  let res = poly_multiplier(&a, &b);

  assert_eq!(
    res,
    TorusPolynomial {
      n: 5,
      coefs: vec![10, 40, 100, 120, 90]
    }
  );
}
