use criterion::{black_box, criterion_group, criterion_main, Criterion};

use tfhe::bootstrapping::{
  boots_sym_decrypt, boots_sym_encrypt, new_default_gate_bootstrapping_parameters,
  new_random_gate_bootstrapping_secret_keyset,
};

/// Benchmarks the encryption function
fn criterion_benchmark(c: &mut Criterion) {
  let message = true;
  let security = 128;
  let params = new_default_gate_bootstrapping_parameters(security);
  let secret_key = new_random_gate_bootstrapping_secret_keyset(&params);
  let encrypted = boots_sym_encrypt(message, &secret_key);
  let mut group = c.benchmark_group("Encryption and decryption");

  group.bench_function("encrypt bit", |b| {
    b.iter_with_large_drop(|| boots_sym_encrypt(black_box(message), &secret_key))
  });

  group.bench_function("decrypt bit", |b| {
    b.iter(|| boots_sym_decrypt(black_box(&encrypted), &secret_key))
  });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
