#include "args.hpp"
#include "nthits.hpp"
#include "user_interface.hpp"
#include "utils.hpp"

#include <argparse/argparse.hpp>
#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <btllib/nthash.hpp>
#include <btllib/seq_reader.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

#define PRINT_EXTRA_BF_STATS                                                                       \
	if (args.verbosity > 1) {                                                                      \
		std::cout << std::endl << "Distinct k-mers Bloom filter stats:" << std::endl;              \
		print_bloom_filter_stats(bf.get_fpr(), args.fpr, bf.get_occupancy());                      \
		std::cout << std::endl << "Intermediate counting Bloom filter stats:" << std::endl;        \
		print_bloom_filter_stats(cbf.get_fpr(), args.fpr, cbf.get_occupancy());                    \
		std::cout << std::endl;                                                                    \
	}

inline unsigned
get_flag(bool long_mode)
{
	if (long_mode) {
		return btllib::SeqReader::Flag::LONG_MODE;
	} else {
		return btllib::SeqReader::Flag::SHORT_MODE;
	}
}

inline std::vector<uint64_t>
load_histogram(const std::string& path)
{
	std::vector<uint64_t> hist;
	std::ifstream hist_file(path);
	std::string freq;
	uint64_t value;
	while (hist_file >> freq >> value) {
		hist.push_back(value);
	}
	return hist;
}

inline void
kmers_cbf(const ProgramArguments& args, size_t bf_size, size_t cbf_size, size_t hits_size)
{
	Timer timer;
	TIMER_START(timer, "Initializing Bloom filters")
	btllib::BloomFilter bf(bf_size, args.num_hashes);
	btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
	btllib::CountingBloomFilter<nthits::cbf_counter_t> hits(hits_size, args.num_hashes);
	TIMER_STOP(timer)
	TIMER_START(timer, "Processing k-mers")
	for (const auto file : args.input_files) {
		btllib::SeqReader reader(file, get_flag(args.long_mode));
#pragma omp parallel shared(reader)
		for (const auto& record : reader) {
			nthits::process(
			    record.seq, args.kmer_length, args.min_count, args.max_count, bf, cbf, hits);
		}
	}
	TIMER_STOP(timer)
	PRINT_EXTRA_BF_STATS
	if (args.verbosity > 0) {
		std::cout << "Output Bloom filter stats:" << std::endl;
		print_bloom_filter_stats(hits.get_fpr(), args.fpr, hits.get_occupancy());
		std::cout << std::endl;
	}
	TIMER_START(timer, "Saving Bloom filter")
	hits.save(args.out_file);
	TIMER_STOP(timer)
}

inline void
seeds_cbf(const ProgramArguments& args, size_t bf_size, size_t cbf_size, size_t hits_size)
{
	Timer timer;
	TIMER_START(timer, "Initializing Bloom filters")
	btllib::BloomFilter bf(bf_size, args.num_hashes);
	btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
	btllib::CountingBloomFilter<nthits::cbf_counter_t> hits(hits_size, args.num_hashes);
	TIMER_STOP(timer)
	TIMER_START(timer, "Processing k-mers")
	for (const auto file : args.input_files) {
		btllib::SeqReader reader(file, get_flag(args.long_mode));
#pragma omp parallel shared(reader)
		for (const auto& record : reader) {
			for (const auto& seed : args.seeds) {
				nthits::process(record.seq, seed, args.min_count, args.max_count, bf, cbf, hits);
			}
		}
	}
	TIMER_STOP(timer)
	PRINT_EXTRA_BF_STATS
	if (args.verbosity > 0) {
		std::cout << "Output Bloom filter stats:" << std::endl;
		print_bloom_filter_stats(hits.get_fpr(), args.fpr, hits.get_occupancy());
		std::cout << std::endl;
	}
	TIMER_START(timer, "Saving Bloom filter")
	hits.save(args.out_file);
	TIMER_STOP(timer)
}

inline void
kmers_table(const ProgramArguments& args, size_t bf_size, size_t cbf_size, size_t hits_size)
{
	Timer timer;
	TIMER_START(timer, "Initializing Bloom filters")
	btllib::BloomFilter bf(bf_size, args.num_hashes);
	btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
	nthits::HitTable hits(hits_size);
	TIMER_STOP(timer)
	TIMER_START(timer, "Processing k-mers")
	for (const auto file : args.input_files) {
		btllib::SeqReader reader(file, get_flag(args.long_mode));
#pragma omp parallel shared(reader)
		for (const auto& record : reader) {
			for (const auto& seed : args.seeds) {
				nthits::process(record.seq, seed, args.min_count, args.max_count, bf, cbf, hits);
			}
		}
	}
	TIMER_STOP(timer)
	PRINT_EXTRA_BF_STATS
	TIMER_START(timer, "Saving Bloom filter")
	hits.save(args.out_file, args.min_count);
	TIMER_STOP(timer)
}

inline void
seeds_table(const ProgramArguments& args, size_t bf_size, size_t cbf_size, size_t hits_size)
{
	Timer timer;
	TIMER_START(timer, "Initializing Bloom filters")
	btllib::BloomFilter bf(bf_size, args.num_hashes);
	btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
	nthits::HitTable hits(hits_size);
	TIMER_STOP(timer)
	TIMER_START(timer, "Processing k-mers")
	for (const auto file : args.input_files) {
		btllib::SeqReader reader(file, get_flag(args.long_mode));
#pragma omp parallel shared(reader)
		for (const auto& record : reader) {
			nthits::process(
			    record.seq, args.kmer_length, args.min_count, args.max_count, bf, cbf, hits);
		}
	}
	TIMER_STOP(timer)
	PRINT_EXTRA_BF_STATS
	TIMER_START(timer, "Saving Bloom filter")
	hits.save(args.out_file, args.min_count);
	TIMER_STOP(timer)
}

int
main(int argc, char** argv)
{
	print_logo();

	auto args = ProgramArguments(argc, argv);

	if (args.verbosity > 1) {
		args.print();
	}

	omp_set_num_threads(args.num_threads);

	auto hist = load_histogram(args.histogram_path);
	size_t hit_count;
	unsigned given_hit_cap = args.min_count;
	nthits::get_thresholds(hist, args.solid, hit_count, args.min_count, args.max_count);
	bool hit_cap_changed = given_hit_cap != args.min_count;

	size_t bf_size, cbf_size, hit_size;
	bf_size = hist[1] * 7 / 8;
	cbf_size = (hist[1] - hist[2]) * 6;
	if (args.out_is_filter()) {
		hit_size = nthits::get_bf_size(hit_count, args.num_hashes, args.seeds.size(), args.fpr);
	} else {
		hit_size = hit_count * 3;
	}
	if (args.verbosity > 0) {
		print_updated_params(
		    hit_count,
		    hist[1],
		    bf_size,
		    args.min_count,
		    hit_cap_changed,
		    cbf_size,
		    args.out_is_filter(),
		    hit_size,
		    args.verbosity);
	}

	if (args.out_type == OutputType::COUNTING_BLOOM_FILTER && !args.using_seeds()) {
		kmers_cbf(args, bf_size, cbf_size, hit_size);
	} else if (args.out_type == OutputType::COUNTING_BLOOM_FILTER && args.using_seeds()) {
		seeds_cbf(args, bf_size, cbf_size, hit_size);
	} else if (args.out_type == OutputType::HIT_TABLE && !args.using_seeds()) {
		kmers_table(args, bf_size, cbf_size, hit_size);
	} else if (args.out_type == OutputType::HIT_TABLE && args.using_seeds()) {
		seeds_table(args, bf_size, cbf_size, hit_size);
	}

	return 0;
}
