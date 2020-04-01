use crate::polynomial::{IntPolynomial, Polynomial, TorusPolynomial};

/// Idea:
/// we may want to represent an element x of the real torus by
/// the integer rint(2^32.x) modulo 2^32
///  -- addition, subtraction and integer combinations are native operation
///  -- modulo 1 is mapped to mod 2^32, which is also native!
/// This looks much better than using float/doubles, where modulo 1 is not
/// natural at all.
pub type Torus32 = i32;

/// 2 ^ 32
const TWO_32: f64 = 4_294_967_296f64;

pub(crate) trait Modulo<RHS = Self> {
  type Output;
  fn modulo(&self, rhs: RHS) -> Self::Output;
}

impl<T: num_traits::PrimInt + num_traits::Signed> Modulo<T> for T {
  type Output = T;
  fn modulo(&self, rhs: T) -> Self::Output {
    let r = *self % rhs;
    if r < T::zero() {
      r + rhs.abs()
    } else {
      r
    }
  }
}

/// Gaussian sample centered in message, with standard deviation sigma
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
  (d * TWO_32).trunc() as Torus32
}

pub(crate) fn torus_32_to_f64(x: Torus32) -> f64 {
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
pub(crate) const fn mod_switch_to_torus32(mu: i32, message_size: i32) -> Torus32 {
  // Width of each interval
  let interval: u64 = ((1u64 << 63) / message_size as u64) * 2;

  // Overflowed here, using wrapping mul to imitate C++ behavior
  let phase64: u64 = (mu as u64).wrapping_mul(interval);

  // Floor to the nearest multiples of interval
  (phase64 >> 32) as Torus32
}

// Used to approximate the phase to the nearest message possible in the message space
// The constant Msize will indicate on which message space we are working (how many messages possible)
//
// "travailler sur 63 bits au lieu de 64, car dans nos cas pratiques, c'est plus précis"
pub(crate) fn mod_switch_from_torus32(phase: Torus32, message_size: i32) -> i32 {
  // width of each intervall
  let interval: u64 = ((1u64 << 63) / message_size as u64) * 2;
  let half_interval: u64 = interval / 2;
  let phase64: u64 = ((phase as u64) << 32).wrapping_add(half_interval);
  //floor to the nearest multiples of interv
  (phase64 / interval) as i32
}

/// Returns the phase of the message in the half-open range `[- message_space / 2, message_space / 2)`
/// in a space with `message_space` possible elements.
/// # Note
/// The `message_space` must be a power of 2.
pub const fn encode_message(mu: i32, message_space: i32) -> Torus32 {
  let log2 = message_space.trailing_zeros();
  mu << (32 - log2)
}

/// Approximates the phase to the nearest possible message in the message space `[- message_space / 2, message_space / 2)`
/// in a space with `message_space` possible elements.
/// # Panics
/// Panics if the `message_space` isn't a power of 2.
pub(crate) fn decode_message(phase: Torus32, message_space: i32) -> i32 {
  assert!(is_power_of_2(message_space as usize));
  let log2 = message_space.trailing_zeros();
  (phase.wrapping_add(1 << (32 - log2 - 1))) >> (32 - log2)
}

/// Pretty self-descriptive
#[inline]
fn is_power_of_2(x: usize) -> bool {
  (x & (x - 1)) == 0
}

pub(crate) fn torus_polynomial_mul_r(
  result: &mut TorusPolynomial,
  poly1: &IntPolynomial,
  poly2: &TorusPolynomial,
) {
  let res = poly_multiplier(poly1, poly2);

  result.coefs = result
    .coefs
    .iter()
    .zip(res.coefs.iter())
    .map(|(a, b)| a + b)
    .collect();
}

/// Multiplies two polynomials
///
/// ## Warning
/// Inefficient -> `O(n<sup>2</sup>)`
#[cfg(not(feature = "fft"))]
pub(crate) fn poly_multiplier<P1, P2>(a: &P1, b: &P2) -> TorusPolynomial
where
  P1: Polynomial<i32>,
  P2: Polynomial<i32>,
{
  let (a_coefs, b_coefs) = crate::polynomial::match_and_pad(a.coefs().to_vec(), b.coefs().to_vec());
  let degree = a.len() + b.len() - 2;
  let mut coefs = vec![0; a_coefs.len() + b_coefs.len() + 1];

  for i in 0..a.coefs().len() {
    for j in 0..b.coefs().len() {
      coefs[i + j] = coefs[i + j] + a_coefs[i] * b_coefs[j];
    }
  }

  TorusPolynomial::from(
    coefs
      .into_iter()
      .skip(a_coefs.len() - degree)
      .collect::<Vec<i32>>(),
  )
}

/// Multiplies two polynomials using FFT
///
/// Efficient -> `O(n log n)`
#[cfg(feature = "fft")]
#[allow(clippy::many_single_char_names)]
pub(crate) fn poly_multiplier<P1, P2>(a: &P1, b: &P2) -> TorusPolynomial
where
  P1: Polynomial<i32>,
  P2: Polynomial<i32>,
{
  // Algorithm found at https://math.stackexchange.com/questions/764727/concrete-fft-polynomial-multiplication-example
  use num_traits::Zero;
  use rustfft::num_complex::Complex;
  use rustfft::FFTplanner;

  let degree = a.len() + b.len() - 2;
  let mut p: Vec<_> = a
    .iter()
    .rev()
    .map(|x| f64::from(*x))
    .chain(std::iter::repeat(0f64).take(b.len() - 1))
    .map(|x| Complex::new(x, 0f64))
    .collect();

  let power = p.len().next_power_of_two();
  if power != p.len() {
    // Extend the polynomial to a power of 2 length
    p.extend(
      std::iter::repeat(Complex::<f64>::zero())
        .take(power - p.len())
        .collect::<Vec<Complex<f64>>>(),
    );
  }

  let mut q: Vec<Complex<f64>> = b
    .iter()
    .rev()
    .map(|x| f64::from(*x))
    .chain(std::iter::repeat(0f64).take(a.len() - 1))
    .map(|x| Complex::new(x, 0f64))
    .collect();

  let power = q.len().next_power_of_two();
  if power != q.len() {
    // Extend the polynomial to a power of 2 length
    q.extend(
      std::iter::repeat(num_traits::identities::Zero::zero())
        .take(power - q.len())
        .collect::<Vec<Complex<f64>>>(),
    );
  }

  let mut p_out = vec![num_traits::identities::Zero::zero(); p.len()];
  let mut q_out = vec![num_traits::identities::Zero::zero(); q.len()];

  // Create a FFT planner for a FFT
  let mut planner = FFTplanner::new(false);
  let fft = planner.plan_fft(p.len());
  fft.process(&mut p, &mut p_out);
  fft.process(&mut q, &mut q_out);

  let mut r: Vec<_> = p_out
    .into_iter()
    .zip(q_out.into_iter())
    .map(|(p_c, q_c)| (p_c / (p.len() as f64).sqrt()) * (q_c / (q.len() as f64).sqrt()))
    .collect();

  let mut res = vec![num_traits::identities::Zero::zero(); r.len()];

  // Create a FFT planner for the inverse FFT
  let fft = FFTplanner::new(true).plan_fft(r.len());
  fft.process(&mut r, &mut res);

  let coefs: Vec<i32> = res
    .into_iter()
    .take(degree + 1)
    .map(|x| x.re.round() as i32)
    .rev()
    .collect();

  TorusPolynomial::from(&coefs)
}

/// X^{a} * source
pub(crate) fn torus_polynomial_mul_by_xai(a: i32, source: &TorusPolynomial) -> TorusPolynomial {
  let n = source.coefs.len() as i32;
  assert!(a >= 0 && a < 2 * n, "a: {}, n: {}, n * 2: {}", a, n, n * 2);
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

  TorusPolynomial::from(coefs)
}

// result = (X^{a}-1)*source
pub(crate) fn torus_polynomial_mul_by_xai_minus_one(
  a: i32,
  source: &TorusPolynomial,
) -> TorusPolynomial {
  let n = source.coefs.len() as i32;
  let mut coefs = vec![0; source.coefs.len()];
  assert!(a >= 0 && a < 2 * n, "{} >= 0 && {} < {}", a, a, 2 * n);

  if a < n {
    for i in 0..a {
      // So that i-a<0 (French: Sur que ...)
      // Overflowed here, using wrapping sub to imitate C++ behavior
      coefs[i as usize] =
        (-source.coefs[(i - a + n) as usize]).wrapping_sub(source.coefs[i as usize]);
    }
    for i in a..n {
      // So that N>i-a>=0 (French: Sur que ...)
      // Overflowed here, using wrapping sub to imitate C++ behavior
      coefs[i as usize] = source.coefs[(i - a) as usize].wrapping_sub(source.coefs[i as usize]);
    }
  } else {
    let aa = a - n;
    for i in 0..aa {
      // So that i-a<0 (French: Sur que ...)
      // Overflowed here, using wrapping sub to imitate C++ behavior
      coefs[i as usize] =
        source.coefs[(i - aa + n) as usize].wrapping_sub(source.coefs[i as usize]);
    }
    for i in aa..n {
      // So that N>i-a>=0 (French: Sur que ...)
      // Overflowed here, using wrapping sub to imitate C++ behavior
      coefs[i as usize] = (-source.coefs[(i - aa) as usize]).wrapping_sub(source.coefs[i as usize]);
    }
  }

  TorusPolynomial::from(coefs)
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::polynomial::Polynomial;
  use rand::distributions::Distribution;

  /// This function return the absolute value of the (centered) fractional part of `d`
  /// i.e. the distance between `d` and its nearest integer
  pub(crate) fn abs_frac(d: f64) -> f64 {
    (d - d.round()).abs()
  }

  /// Unsure what this does, but it works
  fn anticyclic_get(tab: &[i32], a: i32, n: i32) -> i32 {
    let agood = ((a % (2 * n)) + (2 * n)) % (2 * n);
    if agood < n {
      tab[agood as usize]
    } else {
      -tab[(agood - n) as usize]
    }
  }

  #[test]
  fn test_torus_polynomial_mul_by_xai_minus_one() {
    const NB_TRIALS: i32 = 50;
    const DIMENSIONS: [i32; 4] = [500, 750, 1024, 2000];

    let mut rng = rand::thread_rng();
    for &n in DIMENSIONS.iter() {
      for trial in 0..NB_TRIALS {
        let d = rand_distr::Uniform::new(i32::min_value(), i32::max_value());
        let a = (d.sample(&mut rng) % 1_000_000) - 500_000;
        let ai = ((a % (2 * n)) + (2 * n)) % (2 * n);

        // Fill the polynomial with random coefs
        let pola = TorusPolynomial::uniform(n as usize);
        let polb = torus_polynomial_mul_by_xai_minus_one(ai, &pola);

        for j in 0..n {
          assert_eq!(
            polb.coefs[j as usize],
            // Overflowed here, using wrapping sub to imitate C++ behavior
            anticyclic_get(&pola.coefs, j - ai, n).wrapping_sub(anticyclic_get(&pola.coefs, j, n))
          );
        }
      }
    }
  }

  #[test]
  fn test_torus_polynomial_mul_by_xai() {
    const NB_TRIALS: i32 = 50;
    const DIMENSIONS: [i32; 4] = [500, 750, 1024, 2000];

    let mut rng = rand::thread_rng();
    for &n in DIMENSIONS.iter() {
      for trial in 0..NB_TRIALS {
        let d = rand_distr::Uniform::new(i32::min_value(), i32::max_value());
        let a = (d.sample(&mut rng) % 1_000_000) - 500_000;
        let ai = ((a % (2 * n)) + (2 * n)) % (2 * n);

        // Fill the polynomial with random coefs
        let pola = TorusPolynomial::uniform(n as usize);
        let polb = torus_polynomial_mul_by_xai(ai, &pola);

        for j in 0..n {
          assert_eq!(
            polb.coefs[j as usize],
            anticyclic_get(&pola.coefs, j - ai, n)
          );
        }
      }
    }
  }

  #[test]
  fn test_float_to_torus_32() {
    assert_eq!(0, f64_to_torus_32(0f64));
    assert_eq!(1 << 31, f64_to_torus_32(0.5));
    assert_eq!(1 << 31, f64_to_torus_32(-0.5));
    assert_eq!(1 << 30, f64_to_torus_32(0.25));
    assert_eq!(0xC0000000, f64_to_torus_32(-0.25) as u32);
  }

  #[test]
  fn test_approximate_phase() {
    let mut rng = rand::thread_rng();
    let d = rand_distr::Uniform::new(i32::min_value(), i32::max_value());
    for i in 2..200 {
      let v = d.sample(&mut rng);
      let w = approximate_phase(v, i);
      let dv = torus_32_to_f64(v);
      let dw = torus_32_to_f64(w);
      assert!(
        abs_frac(dv - dw) <= 1f64 / (2f64 * i as f64) + 1e-40,
        "{} <= {}",
        abs_frac(dv - dw),
        1f64 / (2f64 * i as f64) + 1e-40
      );
      assert!(
        abs_frac(i as f64 * dw) <= i as f64 * 1e-9,
        "{} <= {}",
        abs_frac(i as f64 * dw),
        i as f64 * 1e-9
      );
    }
  }

  /// Modular gaussian distribution of standard deviation sigma centered on the message, on the Torus32
  #[test]
  fn test_gaussian_torus_32() {
    const MESSAGE1: Torus32 = 123456789;
    const MESSAGE2: Torus32 = 987654321;
    let reps1 = gaussian32(MESSAGE1, 0f64);
    let reps2 = gaussian32(MESSAGE2, 0f64);
    assert_eq!(MESSAGE1, reps1);
    assert_eq!(MESSAGE2, reps2);
    let reps1 = gaussian32(MESSAGE1, 0.01);
    let reps2 = gaussian32(MESSAGE2, 0.5);
    assert_ne!(MESSAGE1, reps1);
    assert_ne!(MESSAGE2, reps2);
    assert!(
      (MESSAGE1 - reps1).abs() <= 80_000_000,
      "{} <= {}",
      (MESSAGE1 - reps1).abs(),
      80_000_000
    );
  }

  #[test]
  /// Converts mu/Msize to a Torus32 for mu in [0,Msize[
  fn test_mod_switch_to_torus_32() {
    let mut rng = rand::thread_rng();
    let d = rand_distr::Uniform::new(i32::min_value(), i32::max_value());
    for i in 2..200 {
      let j = d.sample(&mut rng) % i;
      let v: Torus32 = mod_switch_to_torus32(j, i);

      let dv = torus_32_to_f64(v);
      assert!(
        abs_frac(dv - j as f64 / i as f64) <= 1f64 / (2f64 * i as f64) + 1e-40,
        "{} <= {}",
        abs_frac(dv - j as f64 / i as f64),
        1f64 / (2f64 * i as f64) + 1e-40
      );
    }
  }

  #[test]
  fn test_mod_switch_from_torus_32() {
    let mut rng = rand::thread_rng();
    let d = rand_distr::Uniform::new(i32::min_value(), i32::max_value());
    for i in 2..200 {
      let v: Torus32 = d.sample(&mut rng);
      let w = mod_switch_from_torus32(v, i);
      let dv = torus_32_to_f64(v);
      let dw = w as f64 / i as f64;
      assert!(
        abs_frac(dv - dw) <= 1f64 / (2f64 * i as f64) + 1e-40,
        "{} <= {}",
        abs_frac(dv - dw),
        1f64 / (2f64 * i as f64) + 1e-40
      );
    }
  }

  #[test]
  fn test_poly_multiplier() {
    let a = IntPolynomial::from(vec![10, 20, 30]);
    let b = IntPolynomial::from(vec![1, 2, 3]);

    let res = poly_multiplier(&a, &b);
    assert_eq!(res, TorusPolynomial::from(vec![10, 40, 100, 120, 90]));
  }
}
