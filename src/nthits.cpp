#include "nthits.hpp"
#include "utils.hpp"

#include <btllib/nthash.hpp>
#include <fstream>

namespace {
#define PROCESS_HIT_TABLE                                                                          \
	while (nth.roll()) {                                                                           \
		if (!distincts.contains_insert(nth.hashes())) {                                            \
			std::string current_kmer = seq.substr(nth.get_pos(), distincts.get_k());               \
			to_canonical(current_kmer);                                                            \
			if (hit_cap > 1) {                                                                     \
				if (cbf.insert_thresh_contains(nth.hashes(), hit_cap - 1)) {                       \
					hit_table.insert(nth.hashes()[0], current_kmer);                               \
				}                                                                                  \
			} else {                                                                               \
				hit_table.insert(nth.hashes()[0], current_kmer);                                   \
			}                                                                                      \
		}                                                                                          \
	}

#define PROCESS_HIT_FILTER                                                                         \
	while (nth.roll()) {                                                                           \
		if (distincts.contains_insert(nth.hashes())) {                                             \
			if (hit_cap > 1) {                                                                     \
				if (cbf.insert_thresh_contains(nth.hashes(), hit_cap - 1))                         \
					hit_filter.insert(nth.hashes());                                               \
			} else {                                                                               \
				hit_filter.insert(nth.hashes());                                                   \
			}                                                                                      \
		}                                                                                          \
	}
}

void
nthits::HitTable::insert(const uint64_t hash_value, const std::string& kmer)
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

void
nthits::HitTable::save(const std::string& file_path, const unsigned hit_cap)
{
	std::ofstream outFile(file_path);
	for (size_t i = 0; i < table_size; i++)
		if (entries[i].count != 0)
			outFile << entries[i].kmer << "\t" << entries[i].count + hit_cap << std::endl;
	outFile.close();
}

void
nthits::process(
    const std::string& seq,
    const unsigned hit_cap,
    btllib::KmerBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    HitTable& hit_table)
{
	btllib::NtHash nth(seq, cbf.get_hash_num(), distincts.get_k());
	PROCESS_HIT_TABLE
}

void
nthits::process(
    const std::string& seq,
    const unsigned hit_cap,
    btllib::SeedBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    HitTable& hit_table)
{
	btllib::SeedNtHash nth(seq, distincts.get_seeds(), cbf.get_hash_num(), distincts.get_k());
	PROCESS_HIT_TABLE
}

void
nthits::process(
    const std::string& seq,
    const unsigned hit_cap,
    btllib::KmerBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    btllib::KmerBloomFilter& hit_filter)
{
	btllib::NtHash nth(seq, cbf.get_hash_num(), distincts.get_k());
	PROCESS_HIT_FILTER
}

void
nthits::process(
    const std::string& seq,
    const unsigned hit_cap,
    btllib::SeedBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    btllib::SeedBloomFilter& hit_filter)
{
	btllib::SeedNtHash nth(seq, distincts.get_seeds(), cbf.get_hash_num(), distincts.get_k());
	PROCESS_HIT_FILTER
}
