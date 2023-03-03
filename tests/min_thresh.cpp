#define CATCH_CONFIG_MAIN

#include <btllib/seq_reader.hpp>
#include <btllib/util.hpp>
#include <catch2/catch_test_macros.hpp>
#include <iostream>

#include "nthits_min.hpp"
#include "tests.hpp"

TEST_CASE("BF, cmin=2", "[nthits]")
{
  auto data = load_test_data();
  Validator v(data, 2);
  btllib::BloomFilter bf(BF_SIZE, NUM_HASHES);
  btllib::CountingBloomFilter8 cbf(BF_SIZE, NUM_HASHES);
  btllib::KmerBloomFilter hits(BF_SIZE, NUM_HASHES, v.get_kmer_length());
  for (const auto& seq : data) {
    nthits::find_hits(seq, v.get_kmer_length(), 2, bf, cbf, hits);
  }
  REQUIRE(v.check_kmer_bf(hits));
}