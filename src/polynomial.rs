use crate::numerics::Torus32;
use rand::distributions::Distribution;

#[derive(Clone, Debug, PartialEq)]
pub struct IntPolynomial {
  pub(crate) coefs: Vec<i32>,
}

impl IntPolynomial {
  pub(crate) fn new(n: i32) -> Self {
    Self {
      coefs: vec![0; n as usize],
    }
  }

  /// Euclidean norm, squared
  pub(crate) fn norm_squared(&self) -> f64 {
    self.coefs.iter().map(|c| (c * c) as f64).sum::<f64>()
  }
}

impl From<&[i32]> for IntPolynomial {
  fn from(s: &[i32]) -> Self {
    Self { coefs: s.to_vec() }
  }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct TorusPolynomial {
  pub(crate) coefs: Vec<Torus32>,
}

impl TorusPolynomial {
  pub(crate) fn new(n: i32) -> Self {
    Self {
      coefs: vec![0; n as usize],
    }
  }

  pub(crate) fn torus_polynomial_uniform(n: i32) -> Self {
    let d = rand_distr::Uniform::new(i32::min_value(), i32::max_value());
    let mut rng = rand::thread_rng();

    Self {
      coefs: (0..n as i32).map(|_| d.sample(&mut rng)).collect(),
    }
  }
}

impl From<&[i32]> for TorusPolynomial {
  fn from(s: &[i32]) -> Self {
    Self { coefs: s.to_vec() }
  }
}

impl std::ops::Add<TorusPolynomial> for TorusPolynomial {
  type Output = Self;
  fn add(self, p: Self) -> Self {
    assert_eq!(self.coefs.len(), p.coefs.len());
    Self {
      coefs: self
        .coefs
        .iter()
        .zip(p.coefs.iter())
        .map(|(a, b)| a + b)
        .collect(),
    }
  }
}
