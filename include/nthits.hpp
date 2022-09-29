#ifndef NTHITS_HPP
#define NTHITS_HPP

#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <omp.h>
#include <stdint.h>
#include <string>

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

	void insert(const uint64_t hash_value, const std::string& kmer);

	void save(const std::string& file_path, const unsigned hit_cap);
};

void
process(
    const std::string& seq,
	const unsigned hit_cap,
    btllib::KmerBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    HitTable& hit_table);

void
process(
    const std::string& seq,
	const unsigned hit_cap,
    btllib::SeedBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    HitTable& hit_table);

void
process(
    const std::string& seq,
	const unsigned hit_cap,
    btllib::KmerBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    btllib::KmerBloomFilter& hits);

void
process(
    const std::string& seq,
	const unsigned hit_cap,
    btllib::SeedBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    btllib::SeedBloomFilter& hits);

}

#endif