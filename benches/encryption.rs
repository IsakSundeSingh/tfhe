use criterion::{black_box, criterion_group, criterion_main, Criterion, Throughput};

use tfhe::encryption::{decrypt, encrypt, generate_keys, Parameters, SecurityLevel};

/// Benchmarks the encryption function
fn criterion_benchmark(c: &mut Criterion) {
  let message = true;
  let params = Parameters::with(SecurityLevel::Bit80);
  let (secret_key, _cloud_key) = generate_keys(&params);

  let encrypted = encrypt(message, &secret_key);
  let mut group = c.benchmark_group("Encryption and decryption");
  group.throughput(Throughput::Bytes(1));
  group.bench_function("encrypt bit", |b| {
    b.iter_with_large_drop(|| encrypt(black_box(message), &secret_key))
  });

  group.bench_function("decrypt bit", |b| {
    b.iter(|| decrypt(black_box(&encrypted), &secret_key))
  });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
