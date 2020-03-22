use criterion::{criterion_group, criterion_main, Criterion};

use tfhe::bootstrapping::{
  new_default_gate_bootstrapping_parameters, new_random_gate_bootstrapping_secret_keyset,
};

/// Benchmarks the key generation with the default parameters
fn key_generation_benchmark(c: &mut Criterion) {
  let security = 128;
  let params = new_default_gate_bootstrapping_parameters(security);
  c.bench_function("Generate key", |b| {
    b.iter(|| new_random_gate_bootstrapping_secret_keyset(&params))
  });
}

criterion_group!(
  name = benches;
  config = Criterion::default().sample_size(10);
  targets = key_generation_benchmark
);
criterion_main!(benches);
