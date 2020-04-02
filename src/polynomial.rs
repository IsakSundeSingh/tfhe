use crate::numerics::{Modulo, Torus32};
use num_traits::{int::PrimInt, Zero};
use rand::distributions::Distribution;

/// Polynomials modulo `X^N + 1` or `X^N - 1`, where `N-1` is the polynomial degree
#[derive(Copy, Clone, Debug, PartialEq)]
pub(crate) enum Cyclicity {
  /// A constant denoting negacyclic polynomial modulus (`X^N + 1`), to be supplied to the polynomial constructor.
  Negacyclic,
  /// A constant denoting cyclic polynomial modulus (`X^N - 1`), to be supplied to the polynomial constructor.
  Cyclic,
}

pub(crate) trait Polynomial<T>:
  std::ops::Add<Self> + std::ops::Mul<Self> + std::ops::Index<usize>
where
  Self: Sized,
{
  /// Returns the coefficients of the polynomial.
  ///
  /// # Note
  /// Polynomials are stored in big endian order, meaning the
  /// coefficient of the highest degree term is first in the given slice,
  /// and the coefficient of the term with degree 0 is the last element.
  fn coefs(&self) -> &[T];

  fn iter(&self) -> std::slice::Iter<T>;

  /// Determines the polynomial modulus' cyclicity.
  fn cyclicity(&self) -> Cyclicity;

  /// Generates a random-generated polynomial of length `n` with uniform distribution of elements
  fn uniform(n: usize) -> Self;

  /// Initialize a polynomial with given coefficients and cyclitiy.
  fn with(coefs: &[T], cyclicity: Cyclicity) -> Self;

  /// Initialize a polynomial of size `n` where all values are zero.
  fn zero(n: usize) -> Self;

  /// Determines the number of elements in the polynomial. Should be equal to its degree.
  fn len(&self) -> usize {
    self.coefs().len()
  }

  /// Determines the degree of the polynomial. Should equal the number of elements.
  fn degree(&self) -> usize {
    self.len()
  }

  /// Euclidean norm, squared
  fn norm_squared(&self) -> f64
  where
    f64: From<<T as std::ops::Mul>::Output>,
    T: PrimInt,
  {
    self.iter().map(|&c| f64::from(c * c)).sum::<f64>()
  }

  /// Determines whether a polynomial is zero.
  fn is_zero(&self) -> bool
  where
    T: Zero + PartialEq,
  {
    self.iter().all(|c| *c == Zero::zero())
  }
}

/// Simple function for ensuring two vectors are equal length.
/// Pads them leftwise with `0` until they are equal.
pub(crate) fn match_and_pad<T>(
  a: Vec<T>,
  b: Vec<T>,
) -> (std::collections::VecDeque<T>, std::collections::VecDeque<T>)
where
  T: Default + Clone,
{
  let mut diff = a.len() as i32 - b.len() as i32;

  let mut a_copy = std::collections::VecDeque::from(a);
  let mut b_copy = std::collections::VecDeque::from(b);

  let default = T::default();

  while diff > 0 {
    b_copy.push_front(default.clone());
    diff -= 1;
  }

  while diff < 0 {
    a_copy.push_front(default.clone());
    diff += 1;
  }
  (a_copy, b_copy)
}

#[derive(Clone, Debug, PartialEq)]
pub struct IntPolynomial {
  pub(crate) coefs: Vec<i32>,
  pub(crate) cyclicity: Cyclicity,
}

impl IntPolynomial {
  pub(crate) fn new(n: i32) -> Self {
    Self {
      coefs: vec![0; n as usize],
      cyclicity: Cyclicity::Negacyclic,
    }
  }
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct TorusPolynomial {
  pub(crate) coefs: Vec<Torus32>,
  pub(crate) cyclicity: Cyclicity,
}

impl TorusPolynomial {
  pub(crate) fn new(n: i32) -> Self {
    Self {
      coefs: vec![0; n as usize],
      cyclicity: Cyclicity::Negacyclic,
    }
  }
}

macro_rules! impl_add {
  ($name: ident) => {
    impl std::ops::Add<$name> for $name {
      type Output = Self;
      fn add(self, p: Self) -> Self {
        #[cfg(debug)]
        {
          if self.len() != p.len() {
            println!(
              "Adding polynomials of differing lengths!: len({}) - len({})",
              self.len(),
              p.len()
            )
          }
        }

        let (self_copy, p_copy) = match_and_pad(self.coefs, p.coefs);

        Self::from(
          self_copy
            .iter()
            .zip(p_copy.iter())
            .map(|(a, b)| a + b)
            .collect::<Vec<_>>(),
        )
      }
    }
  };
}
macro_rules! impl_sub {
  ($name: ident) => {
    impl std::ops::Sub<$name> for $name {
      type Output = Self;
      fn sub(self, p: Self) -> Self {
        #[cfg(debug)]
        {
          if self.len() != p.len() {
            println!(
              "Subtracting polynomials of differing lengths!: len({}) - len({})",
              self.len(),
              p.len()
            )
          }
        }

        let (self_copy, p_copy) = match_and_pad(self.coefs, p.coefs);
        let vals = self_copy
          .iter()
          .zip(p_copy.iter())
          .map(|(a, b)| a - b)
          .collect::<Vec<_>>();
        Self::from(vals)
      }
    }
  };
}

macro_rules! impl_mul {
  ($name: ident) => {
    impl std::ops::Mul<$name> for $name {
      type Output = Self;
      fn mul(self, p: Self) -> Self {
        crate::numerics::poly_multiplier(&self, &p).into()
      }
    }
  };
}

macro_rules! impl_polynomial {
  ($name:ident, $ty: ty) => {
    impl Polynomial<$ty> for $name {
      fn cyclicity(&self) -> Cyclicity {
        self.cyclicity
      }

      #[inline]
      fn coefs(&self) -> &[$ty] {
        &self.coefs
      }

      fn iter(&self) -> std::slice::Iter<$ty> {
        self.coefs.iter()
      }

      fn uniform(n: usize) -> Self {
        Self::from(uniform(n))
      }

      fn with(coefs: &[$ty], cyclicity: Cyclicity) -> Self {
        Self {
          coefs: coefs.to_vec(),
          cyclicity,
        }
      }

      fn zero(n: usize) -> Self {
        Self::from(vec![0; n])
      }
    }
  };
}

macro_rules! impl_from {
  ($name: ident, $ty: ty) => {
    impl<T> From<T> for $name
    where
      T: AsRef<[$ty]>,
    {
      #[inline]
      fn from(s: T) -> Self {
        let coefs = s.as_ref();
        Self {
          coefs: coefs.to_vec(),
          cyclicity: Cyclicity::Negacyclic,
        }
      }
    }
  };
}

macro_rules! impl_from_poly {
  ($name: ident, $ty: ty) => {
    impl From<$ty> for $name {
      #[inline]
      fn from(p: $ty) -> Self {
        Self {
          coefs: p.coefs,
          cyclicity: p.cyclicity,
        }
      }
    }
  };
}

macro_rules! impl_index {
  ($name: ident, $ty: ty) => {
    impl std::ops::Index<usize> for $name {
      type Output = $ty;
      #[inline]
      fn index(&self, idx: usize) -> &$ty {
        &self.coefs[idx]
      }
    }
  };
}

impl_add!(IntPolynomial);
impl_add!(TorusPolynomial);
impl_sub!(IntPolynomial);
impl_sub!(TorusPolynomial);
impl_mul!(IntPolynomial);
impl_mul!(TorusPolynomial);
impl_from!(IntPolynomial, i32);
impl_from!(TorusPolynomial, i32);
impl_from_poly!(TorusPolynomial, IntPolynomial);
impl_from_poly!(IntPolynomial, TorusPolynomial);
impl_index!(IntPolynomial, i32);
impl_index!(TorusPolynomial, i32);
impl_polynomial!(IntPolynomial, i32);
impl_polynomial!(TorusPolynomial, i32);

/// Simply generates a vector of length `n` with uniform distribution
fn uniform(n: usize) -> Vec<i32> {
  let d = rand_distr::Uniform::new(i32::min_value(), i32::max_value());
  let mut rng = rand::thread_rng();
  (0..n as i32).map(|_| d.sample(&mut rng)).collect()
}

/// Multiply the polynomial by `x^power`.
/// If `power` lies outside `[0, 2 * N)` where `N-1` is the maximum degree of the polynomial,
/// a modulo `2 * N` will be taken.
pub(crate) fn mul_by_monomial(p: IntPolynomial, power: i32) -> IntPolynomial {
  if power == 0 {
    return p;
  }

  let n = p.len() as i32;
  // let newpower = (dbg!(power).modulo(dbg!(n) as i32)) as usize;
  // let power = dbg!(newpower);
  let cycle: bool = (power / n as i32).modulo(2) == 1;
  let power = power % n as i32;
  let mut coefs = vec![0; n as usize];

  let shift_first = matches!(p.cyclicity(), Cyclicity::Negacyclic) && (!cycle);
  let shift_last = matches!(p.cyclicity(), Cyclicity::Negacyclic) && cycle;

  for j in 0..power {
    coefs[(j.modulo(n)) as usize] = if shift_first {
      -p[((n - power + j).modulo(n)) as usize]
    } else {
      p[((n - power + j).modulo(n)) as usize]
    }
  }

  for j in power..n {
    coefs[(j.modulo(n)) as usize] = if shift_last {
      -p[((j - power).modulo(n)) as usize]
    } else {
      p[((j - power).modulo(n)) as usize]
    }
  }
  IntPolynomial::with(&coefs, p.cyclicity())
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_sub_different_sizes() {
    fn t(a: &[Torus32], b: &[Torus32], expected: &[Torus32]) {
      let a = TorusPolynomial::from(a);
      let b = TorusPolynomial::from(b);
      let res = a - b;
      assert_eq!(res, TorusPolynomial::from(expected));
    }

    t(&[1, 2, 3, 4], &[1, 1, 2, 3, 4], &[-1, 0, 0, 0, 0]);
    t(&[1, 1, 2, 3, 4], &[1, 2, 3, 4], &[1, 0, 0, 0, 0]);
    t(&[1, 2, 3, 4], &[1, 2, 3, 4], &[0, 0, 0, 0]);
  }

  #[test]
  fn test_mul_by_monomial() {
    let coefs: Vec<i32> = (0..10).collect();
    let m = 21;
    let mr = m as u8;
    let p = IntPolynomial::from(coefs);
    let s = -1;

    assert_eq!(
      mul_by_monomial(p.clone(), 4),
      IntPolynomial::from((6..10).map(|x| x * s).chain(0..6).collect::<Vec<_>>())
    );

    assert_eq!(
      mul_by_monomial(p.clone(), 10),
      IntPolynomial::from((0..10).map(|x| x * s).collect::<Vec<_>>())
    );

    assert_eq!(
      mul_by_monomial(p.clone(), 14),
      IntPolynomial::from((6..10).chain((0..6).map(|x| x * s)).collect::<Vec<_>>())
    );

    assert_eq!(
      mul_by_monomial(p.clone(), 20),
      IntPolynomial::from((0..10).collect::<Vec<_>>())
    );

    assert_eq!(
      mul_by_monomial(p.clone(), 24),
      IntPolynomial::from(((6..10).map(|x| x * s)).chain(0..6).collect::<Vec<_>>())
    );

    // @test mul_by_monomial(p, -4) == Polynomial(mtp.([4:9; s .* (0:3)]), pm)
    // @test mul_by_monomial(p, -10) == Polynomial(mtp.([s .* (0:9);]), pm)
    // @test mul_by_monomial(p, -14) == Polynomial(mtp.([s .* (4:9); 0:3]), pm)
    // @test mul_by_monomial(p, -20) == Polynomial(mtp.([0:9;]), pm)
    // @test mul_by_monomial(p, -24) == Polynomial(mtp.([4:9; s .* (0:3)]), pm)
  }
}
