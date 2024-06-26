#ifndef NTHITS_MIN_HPP
#define NTHITS_MIN_HPP

#include "hit_table.hpp"
#include "nthits.hpp"

#define PROCESS_MIN(HIT_INIT, HIT_INSERT)                                                          \
  while (nthash.roll()) {                                                                          \
    if (!bf.contains_insert(nthash.hashes())) {                                                    \
      continue;                                                                                    \
    }                                                                                              \
    if (min_count <= 2) {                                                                          \
      HIT_INSERT                                                                                   \
      continue;                                                                                    \
    }                                                                                              \
    const unsigned count = cbf.insert_contains(nthash.hashes()) + 1;                               \
    if (count == min_count) {                                                                      \
      HIT_INIT                                                                                     \
    } else if (count > min_count) {                                                                \
      HIT_INSERT                                                                                   \
    }                                                                                              \
  }

namespace nthits {

inline void
find_hits(const std::string& seq,
          const unsigned kmer_length,
          const unsigned min_count,
          btllib::BloomFilter& bf,
          btllib::CountingBloomFilter<cbf_counter_t>& cbf,
          btllib::KmerBloomFilter& hit_filter)
{
  SEQ_LEN_GUARD(kmer_length)
  btllib::NtHash nthash(seq, bf.get_hash_num(), kmer_length);
  PROCESS_MIN(hit_filter.insert(nthash.hashes());, hit_filter.insert(nthash.hashes());)
}

inline void
find_hits(const std::string& seq,
          const std::string& seed,
          const unsigned min_count,
          btllib::BloomFilter& bf,
          btllib::CountingBloomFilter<cbf_counter_t>& cbf,
          btllib::SeedBloomFilter& hit_filter)
{
  SEQ_LEN_GUARD(seed.size())
  btllib::SeedNtHash nthash(seq, { seed }, bf.get_hash_num(), seed.size());
  PROCESS_MIN(hit_filter.insert(nthash.hashes());, hit_filter.insert(nthash.hashes());)
}

inline void
find_hits(const std::string& seq,
          const unsigned kmer_length,
          const unsigned min_count,
          btllib::BloomFilter& bf,
          btllib::CountingBloomFilter<cbf_counter_t>& cbf,
          btllib::KmerCountingBloomFilter<cbf_counter_t>& hit_filter)
{
  SEQ_LEN_GUARD(kmer_length)
  btllib::NtHash nthash(seq, bf.get_hash_num(), kmer_length);
  PROCESS_MIN(hit_filter.insert(nthash.hashes(), min_count);, hit_filter.insert(nthash.hashes());)
}

inline void
find_hits(const std::string& seq,
          const std::string& seed,
          const unsigned min_count,
          btllib::BloomFilter& bf,
          btllib::CountingBloomFilter<cbf_counter_t>& cbf,
          btllib::CountingBloomFilter<cbf_counter_t>& hit_filter)
{
  SEQ_LEN_GUARD(seed.size())
  btllib::SeedNtHash nthash(seq, { seed }, bf.get_hash_num(), seed.size());
  PROCESS_MIN(hit_filter.insert(nthash.hashes(), min_count);, hit_filter.insert(nthash.hashes());)
}

inline void
find_hits(const std::string& seq,
          const unsigned kmer_length,
          const unsigned min_count,
          btllib::BloomFilter& bf,
          btllib::CountingBloomFilter<cbf_counter_t>& cbf,
          HitTable& hit_table)
{
  SEQ_LEN_GUARD(kmer_length)
  btllib::NtHash nthash(seq, bf.get_hash_num(), kmer_length);
  PROCESS_MIN(
    {
      std::string kmer = seq.substr(nthash.get_pos(), nthash.get_k());
      hit_table.insert(nthash.hashes()[0], kmer);
    },
    {
      std::string kmer = seq.substr(nthash.get_pos(), nthash.get_k());
      hit_table.insert(nthash.hashes()[0], kmer, min_count);
    })
}

}

#endif