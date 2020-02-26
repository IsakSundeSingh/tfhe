use crate::tlwe::Torus32;

// Gaussian sample centered in message, with standard deviation sigma
pub(crate) fn gaussian32(message: Torus32, sigma: f64) -> Torus32 {
  use rand::distributions::Distribution;

  // Attention: all the implementation will use the stdev instead of the gaussian fourier param
  let d = rand::distributions::Normal::new(0f64, sigma);
  let mut rng = rand::thread_rng();
  let error = d.sample(&mut rng);
  message + f64_to_torus_32(error)
}

/// # Warning
/// Weird and lossy conversion
pub(crate) fn f64_to_torus_32(d: f64) -> Torus32 {
  const TWO_32: i64 = 4294967296; // 2 ^ 32
  let inner = d - (d as i64) as f64;
  let x = (inner * TWO_32 as f64) as i64;
  x as Torus32
  // C++ conversion: // return int32_t(int64_t((d - int64_t(d))*_two32));
}
