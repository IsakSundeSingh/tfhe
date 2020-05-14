use criterion::{criterion_group, criterion_main, Criterion};

use tfhe::encryption::{generate_keys, Parameters, SecurityLevel};

/// Benchmarks the key generation with the default parameters
fn key_generation_benchmark(c: &mut Criterion) {
  let params = Parameters::with(SecurityLevel::Bit80);
  c.bench_function("Generate key", |b| b.iter(|| generate_keys(&params)));
}

criterion_group!(
  name = benches;
  config = Criterion::default().sample_size(10);
  targets = key_generation_benchmark
);
criterion_main!(benches);
