#ifndef NTHITS_HPP
#define NTHITS_HPP

#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <btllib/nthash.hpp>
#include <string>

#include "hit_table.hpp"

#define SEQ_LEN_GUARD(K)                                                                           \
  if (seq.size() < K) {                                                                            \
    return;                                                                                        \
  }

#define POPULATE_BF                                                                                \
  while (nthash.roll()) {                                                                          \
    bf.insert(nthash.hashes());                                                                    \
  }

namespace nthits {

const std::string VERSION = "@PROJECT_VERSION@";

using cbf_counter_t = uint8_t;

inline void
find_hits(const std::string& seq, const unsigned kmer_length, btllib::KmerBloomFilter& bf)
{
  SEQ_LEN_GUARD(kmer_length)
  btllib::NtHash nthash(seq, bf.get_hash_num(), kmer_length);
  POPULATE_BF
}

inline void
find_hits(const std::string& seq, const std::string& seed, btllib::BloomFilter& bf)
{
  SEQ_LEN_GUARD(seed.size())
  btllib::SeedNtHash nthash(seq, { seed }, bf.get_hash_num(), seed.size());
  POPULATE_BF
}

inline void
find_hits(const std::string& seq,
          const unsigned kmer_length,
          btllib::KmerCountingBloomFilter<cbf_counter_t>& bf)
{
  SEQ_LEN_GUARD(kmer_length)
  btllib::NtHash nthash(seq, bf.get_hash_num(), kmer_length);
  POPULATE_BF
}

inline void
find_hits(const std::string& seq,
          const std::string& seed,
          btllib::CountingBloomFilter<cbf_counter_t>& bf)
{
  SEQ_LEN_GUARD(seed.size())
  btllib::SeedNtHash nthash(seq, { seed }, bf.get_hash_num(), seed.size());
  POPULATE_BF
}

inline void
find_hits(const std::string& seq, const unsigned kmer_length, HitTable& hits)
{
  SEQ_LEN_GUARD(kmer_length)
  btllib::NtHash nthash(seq, 1, kmer_length);
  while (nthash.roll()) {
    std::string kmer = seq.substr(nthash.get_pos(), nthash.get_k());
    hits.insert(nthash.hashes()[0], kmer);
  }
}

}

#endif