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

#define POPULATE_KMERS(BF, CBF, HITS)                                                              \
	for (const auto& file_path : args.input_files) {                                               \
		btllib::SeqReader reader(file_path, seq_reader_mode);                                      \
		_Pragma("omp parallel shared(reader)") for (const auto record : reader)                    \
		{                                                                                          \
			if (record.seq.size() < args.kmer_length) {                                            \
				continue;                                                                          \
			}                                                                                      \
			nthits::process(                                                                       \
			    record.seq,                                                                        \
			    args.kmer_length,                                                                  \
			    args.num_hashes,                                                                   \
			    args.min_count,                                                                    \
			    args.max_count,                                                                    \
			    BF,                                                                                \
			    CBF,                                                                               \
			    HITS);                                                                             \
		}                                                                                          \
	}

#define POPULATE_SEEDS(BF, CBF, HITS)                                                              \
	for (const auto& file_path : args.input_files) {                                               \
		btllib::SeqReader reader(file_path, seq_reader_mode);                                      \
		_Pragma("omp parallel shared(reader)") for (const auto record : reader)                    \
		{                                                                                          \
			for (const auto& seed : args.seeds) {                                                  \
				if (record.seq.size() < seed.size()) {                                             \
					continue;                                                                      \
				}                                                                                  \
				nthits::process(                                                                   \
				    record.seq,                                                                    \
				    seed,                                                                          \
				    args.num_hashes,                                                               \
				    args.min_count,                                                                \
				    args.max_count,                                                                \
				    BF,                                                                            \
				    CBF,                                                                           \
				    HITS);                                                                         \
			}                                                                                      \
		}                                                                                          \
	}

#define PRINT_EXTRA_BF_STATS                                                                       \
	if (args.verbosity > 1) {                                                                      \
		std::cout << std::endl << "Distinct k-mers Bloom filter stats:" << std::endl;              \
		print_bloom_filter_stats(bf.get_fpr(), args.fpr, bf.get_occupancy());                      \
		std::cout << std::endl << "Intermediate counting Bloom filter stats" << std::endl;         \
		print_bloom_filter_stats(cbf.get_fpr(), args.fpr, cbf.get_occupancy());                    \
		std::cout << std::endl;                                                                    \
	}

int
main(int argc, char** argv)
{
	auto args = parse_arguments(argc, argv);

	if (args.verbosity > 0) {
		print_logo();
	}

	if (args.verbosity > 1) {
		print_args(args);
	}

	unsigned seq_reader_mode;
	if (args.long_mode) {
		seq_reader_mode = btllib::SeqReader::Flag::LONG_MODE;
	} else {
		seq_reader_mode = btllib::SeqReader::Flag::SHORT_MODE;
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
	if (args.out_bloom) {
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
		    args.out_bloom,
		    hit_size,
		    args.verbosity);
	}

	Timer timer;
	if (args.out_bloom && args.seeds.empty()) {
		TIME_EXECUTION(
		    "Initializing Bloom filters", timer, btllib::BloomFilter bf(bf_size, args.num_hashes);
		    btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
		    btllib::KmerCountingBloomFilter<nthits::cbf_counter_t> hits(
		        hit_size, args.num_hashes, args.kmer_length);)
		TIME_EXECUTION("Processing k-mers", timer, POPULATE_KMERS(bf, cbf, hits))
		PRINT_EXTRA_BF_STATS
		if (args.verbosity > 0) {
			std::cout << "Output Bloom filter stats:" << std::endl;
			print_bloom_filter_stats(hits.get_fpr(), args.fpr, hits.get_occupancy());
			std::cout << std::endl;
		}
		TIME_EXECUTION("Saving Bloom filter", timer, hits.save(args.out_file);)
	} else if (args.out_bloom) {
		TIME_EXECUTION(
		    "Initializing Bloom filters", timer, btllib::BloomFilter bf(bf_size, args.num_hashes);
		    btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
		    btllib::KmerCountingBloomFilter<nthits::cbf_counter_t> hits(
		        hit_size, args.num_hashes, args.seeds[0].size());)
		TIME_EXECUTION("Processing spaced seeds", timer, POPULATE_SEEDS(bf, cbf, hits))
		PRINT_EXTRA_BF_STATS
		if (args.verbosity > 0) {
			std::cout << "Output Bloom filter stats:" << std::endl;
			print_bloom_filter_stats(hits.get_fpr(), args.fpr, hits.get_occupancy());
			std::cout << std::endl;
		}
		TIME_EXECUTION("Saving Bloom filter", timer, hits.save(args.out_file);)
		PRINT_EXTRA_BF_STATS
	} else if (args.seeds.empty()) {
		TIME_EXECUTION(
		    "Initializing Bloom filters and hit table",
		    timer,
		    btllib::BloomFilter bf(bf_size, args.num_hashes);
		    btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
		    nthits::HitTable hits(hit_size);)
		TIME_EXECUTION("Processing k-mers", timer, POPULATE_KMERS(bf, cbf, hits))
		PRINT_EXTRA_BF_STATS
		TIME_EXECUTION("Saving hits table", timer, hits.save(args.out_file, args.min_count);)
	} else {
		TIME_EXECUTION(
		    "Initializing Bloom filters and hit table",
		    timer,
		    btllib::BloomFilter bf(bf_size, args.num_hashes);
		    btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
		    nthits::HitTable hits(hit_size);)
		TIME_EXECUTION("Processing spaced seeds", timer, POPULATE_SEEDS(bf, cbf, hits))
		PRINT_EXTRA_BF_STATS
		TIME_EXECUTION("Saving hits table", timer, hits.save(args.out_file, args.min_count);)
	}

	return 0;
}
