use criterion::{criterion_group, criterion_main, BatchSize, Criterion, Throughput};
use tfhe::{
  encryption::{encrypt, generate_keys, Parameters},
  numerics::encode_message,
  tfhe_bootstrap, SecurityLevel,
};

const SAMPLE_SIZE: usize = 10;

fn bootstrapping_benchmark(c: &mut Criterion) {
  let message = true;
  let params = Parameters::with(SecurityLevel::Bit80);
  let (secret_key, cloud_key) = generate_keys(&params);
  let s = encrypt(message, &secret_key);

  let mut group = c.benchmark_group("");
  group.throughput(Throughput::Bytes(2048));
  group.bench_function("bootstrapping", |b| {
    b.iter_batched(
      || s.clone(),
      |data| tfhe_bootstrap(&cloud_key.bk, encode_message(1, 8), data),
      BatchSize::SmallInput,
    )
  });
}

criterion_group!(
  name = benches;
  config = Criterion::default().sample_size(SAMPLE_SIZE).measurement_time(std::time::Duration::from_secs(75));
  targets = bootstrapping_benchmark
);

criterion_main!(benches);
