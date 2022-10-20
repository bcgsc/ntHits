#ifndef NTHITS_HPP
#define NTHITS_HPP

#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <btllib/nthash.hpp>
#include <fstream>
#include <omp.h>
#include <stdint.h>
#include <string>

#include "utils.hpp"

namespace {
#define PROCESS(NTHASH, HIT_INSERT)                                                                \
	NTHASH                                                                                         \
	while (nth.roll()) {                                                                           \
		if (distincts.contains_insert(nth.hashes())) {                                             \
			if (hit_cap > 1) {                                                                     \
				if (cbf.insert_contains(nth.hashes()) > hit_cap - 1) {                             \
					HIT_INSERT                                                                     \
				}                                                                                  \
			} else {                                                                               \
				HIT_INSERT                                                                         \
			}                                                                                      \
		}                                                                                          \
	}
}

namespace nthits {

// Counting Bloom filter bucket type
using cbf_counter_t = uint8_t;

class HitTable
{
  private:
	struct TableEntry
	{
		std::string kmer;
		unsigned count;
	};

	const size_t table_size;
	TableEntry* entries;
	const unsigned lock_size;
	omp_lock_t* locks;

  public:
	explicit HitTable(const size_t table_size, const unsigned lock_size = 65536)
	  : table_size(table_size)
	  , entries(new TableEntry[table_size])
	  , lock_size(lock_size)
	  , locks(new omp_lock_t[lock_size])
	{
		for (size_t i = 0; i < table_size; i++) {
			entries[i].count = 0;
		}
		for (unsigned i = 0; i < lock_size; i++) {
			omp_init_lock(&locks[i]);
		}
	}

	~HitTable()
	{
		for (unsigned i = 0; i < lock_size; i++)
			omp_destroy_lock(&locks[i]);
		delete[] locks;
		delete[] entries;
	}

	void insert(const uint64_t hash_value, const std::string& kmer)
	{
		uint64_t i = 0, j;
		do {
			j = (hash_value + i) % table_size;
			if (entries[j].kmer == kmer) {
#pragma omp atomic
				++entries[j].count;
			}
			++i;
		} while (i != table_size && entries[j].count != 0);
		if (entries[j].count == 0) {
			omp_set_lock(&locks[(uint16_t)j]);
			entries[j].kmer = kmer;
			++entries[j].count;
			omp_unset_lock(&locks[(uint16_t)j]);
		}
	}

	void save(const std::string& file_path, const unsigned hit_cap)
	{
		std::ofstream outFile(file_path);
		for (size_t i = 0; i < table_size; i++)
			if (entries[i].count != 0)
				outFile << entries[i].kmer << "\t" << entries[i].count + hit_cap << std::endl;
		outFile.close();
	}
};

inline void
process(
    const std::string& seq,
    const unsigned hit_cap,
    btllib::KmerBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    HitTable& hit_table)
{
	PROCESS(btllib::NtHash nth(seq, distincts.get_hash_num(), distincts.get_k());
	        , std::string current_kmer = seq.substr(nth.get_pos(), nth.get_k());
	        to_canonical(current_kmer);
	        hit_table.insert(nth.hashes()[0], current_kmer);)
}

inline void
process(
    const std::string& seq,
    const std::string& seed,
    const unsigned hit_cap,
    btllib::SeedBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    HitTable& hit_table)
{
	PROCESS(
	    btllib::SeedNtHash nth(seq, { seed }, distincts.get_hash_num_per_seed(), distincts.get_k());
	    , std::string current_kmer = seq.substr(nth.get_pos(), nth.get_k());
	    to_canonical(current_kmer);
	    hit_table.insert(nth.hashes()[0], current_kmer);)
}

inline void
process(
    const std::string& seq,
    const unsigned hit_cap,
    btllib::KmerBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    btllib::KmerBloomFilter& hit_filter)
{
	PROCESS(btllib::NtHash nth(seq, hit_filter.get_hash_num(), distincts.get_k());
	        , hit_filter.insert(nth.hashes());)
}

inline void
process(
    const std::string& seq,
    const std::string& seed,
    const unsigned hit_cap,
    btllib::SeedBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    btllib::SeedBloomFilter& hit_filter)
{
	PROCESS(
	    btllib::SeedNtHash nth(seq, { seed }, distincts.get_hash_num_per_seed(), distincts.get_k());
	    , hit_filter.insert(nth.hashes());)
}

}

#endif