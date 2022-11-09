#ifndef NTHITS_MIN_MAX_HPP
#define NTHITS_MIN_MAX_HPP

#include "nthits.hpp"

#define PROCESS_MIN_MAX(HIT_INSERT, HIT_REMOVE)                                                    \
  while (nthash.roll()) {                                                                          \
    if (bf.contains_insert(nthash.hashes())) {                                                     \
      auto count = cbf.insert_thresh_contains(nthash.hashes(), max_count) + 1U;                    \
      if (count == max_count + 1) {                                                                \
        HIT_REMOVE                                                                                 \
      } else if (count >= min_count && count <= max_count) {                                       \
        HIT_INSERT                                                                                 \
      }                                                                                            \
    }                                                                                              \
  }

namespace nthits {

inline void
find_hits(const std::string& seq,
          const unsigned kmer_length,
          const unsigned min_count,
          const unsigned max_count,
          btllib::BloomFilter& bf,
          btllib::CountingBloomFilter<cbf_counter_t>& cbf,
          btllib::KmerCountingBloomFilter<cbf_counter_t>& hit_filter)
{
  SEQ_LEN_GUARD(kmer_length)
  btllib::NtHash nthash(seq, bf.get_hash_num(), kmer_length);
  PROCESS_MIN_MAX(hit_filter.insert(nthash.hashes());, hit_filter.clear(nthash.hashes());)
}

inline void
find_hits(const std::string& seq,
          const std::string& seed,
          const unsigned min_count,
          const unsigned max_count,
          btllib::BloomFilter& bf,
          btllib::CountingBloomFilter<cbf_counter_t>& cbf,
          btllib::CountingBloomFilter<cbf_counter_t>& hit_filter)
{
  SEQ_LEN_GUARD(seed.size())
  btllib::SeedNtHash nthash(seq, { seed }, bf.get_hash_num(), seed.size());
  PROCESS_MIN_MAX(hit_filter.insert(nthash.hashes());, hit_filter.clear(nthash.hashes());)
}

inline void
find_hits(const std::string& seq,
          const unsigned kmer_length,
          const unsigned min_count,
          const unsigned max_count,
          btllib::BloomFilter& bf,
          btllib::CountingBloomFilter<cbf_counter_t>& cbf,
          HitTable& hit_table)
{
  SEQ_LEN_GUARD(kmer_length)
  btllib::NtHash nthash(seq, bf.get_hash_num(), kmer_length);
  PROCESS_MIN_MAX(
    {
      std::string kmer = seq.substr(nthash.get_pos(), nthash.get_k());
      hit_table.insert(nthash.hashes()[0], kmer);
    },
    hit_table.remove(nthash.hashes()[0]);)
}

}

#endif