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
#define PROCESS_MIN_MAX(HIT_INSERT, HIT_REMOVE)                                                    \
	while (nthash.roll()) {                                                                        \
		if (bf.contains_insert(nthash.hashes())) {                                                 \
			auto count = cbf.insert_thresh_contains(nthash.hashes(), max_count) + 1U;              \
			if (count == max_count + 1) {                                                          \
				HIT_REMOVE                                                                         \
			} else if (count >= min_count && count <= max_count) {                                 \
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
		std::string canonical = kmer;
		to_canonical(canonical);
		uint64_t i = 0, j;
		do {
			j = (hash_value + i) % table_size;
			if (entries[j].kmer == canonical) {
#pragma omp atomic
				++entries[j].count;
			}
			++i;
		} while (i != table_size && entries[j].count != 0);
		if (entries[j].count == 0) {
			omp_set_lock(&locks[(uint16_t)j]);
			entries[j].kmer = canonical;
			++entries[j].count;
			omp_unset_lock(&locks[(uint16_t)j]);
		}
	}

	void remove(const uint64_t hash_value)
	{ // TODO
	}

	void save(const std::string& file_path, const unsigned min_count)
	{
		std::ofstream outFile(file_path);
		for (size_t i = 0; i < table_size; i++)
			if (entries[i].count != 0)
				outFile << entries[i].kmer << "\t" << entries[i].count + min_count << std::endl;
		outFile.close();
	}
};

inline void
process(
    const std::string& seq,
    const unsigned kmer_length,
    const unsigned min_count,
    const unsigned max_count,
    btllib::BloomFilter& bf,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    HitTable& hit_table)
{
	btllib::NtHash nthash(seq, bf.get_hash_num(), kmer_length);
	PROCESS_MIN_MAX(std::string kmer = seq.substr(nthash.get_pos(), nthash.get_k());
	                hit_table.insert(nthash.hashes()[0], kmer);
	                , hit_table.remove(nthash.hashes()[0]);)
}

inline void
process(
    const std::string& seq,
    const std::string& seed,
    const unsigned min_count,
    const unsigned max_count,
    btllib::BloomFilter& bf,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    HitTable& hit_table)
{
	btllib::SeedNtHash nthash(seq, { seed }, bf.get_hash_num(), seed.size());
	PROCESS_MIN_MAX(std::string kmer = seq.substr(nthash.get_pos(), nthash.get_k());
	                hit_table.insert(nthash.hashes()[0], kmer);
	                , hit_table.remove(nthash.hashes()[0]);)
}

inline void
process(
    const std::string& seq,
    const unsigned kmer_length,
    const unsigned min_count,
    const unsigned max_count,
    btllib::BloomFilter& bf,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    btllib::CountingBloomFilter<cbf_counter_t>& hit_filter)
{
	btllib::NtHash nthash(seq, bf.get_hash_num(), kmer_length);
	PROCESS_MIN_MAX(hit_filter.insert(nthash.hashes());, hit_filter.clear(nthash.hashes());)
}

inline void
process(
    const std::string& seq,
    const std::string& seed,
    const unsigned min_count,
    const unsigned max_count,
    btllib::BloomFilter& bf,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    btllib::CountingBloomFilter<cbf_counter_t>& hit_filter)
{
	btllib::SeedNtHash nthash(seq, { seed }, bf.get_hash_num(), seed.size());
	PROCESS_MIN_MAX(hit_filter.insert(nthash.hashes());, hit_filter.clear(nthash.hashes());)
}

}

#endif