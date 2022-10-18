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

std::vector<uint64_t>
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

#define POPULATE_KMERS(HITS_CONTAINER)                                                             \
	for (const auto& file_path : args.input_files) {                                               \
		btllib::SeqReader reader(file_path, seq_reader_mode);                                      \
		_Pragma("omp parallel shared(reader)") for (const auto record : reader)                    \
		{                                                                                          \
			if (record.seq.size() < args.kmer_length) {                                            \
				continue;                                                                          \
			}                                                                                      \
			nthits::process(record.seq, args.thresh_min, bf, cbf, HITS_CONTAINER);                 \
		}                                                                                          \
	}

#define POPULATE_SEEDS(HITS_CONTAINER)                                                             \
	for (const auto& file_path : args.input_files) {                                               \
		btllib::SeqReader reader(file_path, seq_reader_mode);                                      \
		_Pragma("omp parallel shared(reader)") for (const auto record : reader)                    \
		{                                                                                          \
			for (const auto& seed : args.seeds) {                                                  \
				if (record.seq.size() < args.kmer_length) {                                        \
					continue;                                                                      \
				}                                                                                  \
				nthits::process(record.seq, seed, args.thresh_min, bf, cbf, HITS_CONTAINER);       \
			}                                                                                      \
		}                                                                                          \
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

	size_t dbf_size, cbf_size, hit_size;
	unsigned num_seeds = args.seeds.size();
	auto hist = load_histogram(args.histogram_path);
	size_t hit_count;
	unsigned given_hit_cap = args.thresh_min;
	nthits::get_thresholds(hist, args.solid, hit_count, args.thresh_min);
	bool hit_cap_changed = given_hit_cap != args.thresh_min;
	dbf_size = nthits::calc_bf_size(hist[1], args.num_hashes, num_seeds, args.fpr);
	cbf_size = nthits::calc_bf_size(hist[1] - hist[2], args.num_hashes, num_seeds, args.fpr);
	if (args.out_bloom) {
		hit_size = nthits::calc_bf_size(hit_count, args.num_hashes, num_seeds, args.fpr);
	} else {
		hit_size = hit_count * 3;
	}
	if (args.verbosity > 0) {
		print_updated_params(
		    hit_count, hist[1], args.thresh_min, hit_cap_changed, args.out_bloom, hit_size);
	}

	Timer timer;
	std::cout << "Running ntHits... " << std::flush;
	timer.start();
	btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
	double bf_fpr = 0, bf_occ = 0;
	if (args.out_bloom && args.seeds.empty()) {
		btllib::KmerBloomFilter bf(dbf_size, args.num_hashes, args.kmer_length);
		btllib::KmerBloomFilter hits_filter(hit_size, args.num_hashes, args.kmer_length);
		POPULATE_KMERS(hits_filter)
		hits_filter.save(args.out_file);
		bf_fpr = hits_filter.get_fpr();
		bf_occ = hits_filter.get_occupancy();
	} else if (args.out_bloom) {
		btllib::SeedBloomFilter bf(dbf_size, args.kmer_length, args.seeds, args.num_hashes);
		btllib::SeedBloomFilter hits_filter(
		    hit_size, args.kmer_length, args.seeds, args.num_hashes);
		POPULATE_SEEDS(hits_filter)
		hits_filter.save(args.out_file);
		bf_fpr = hits_filter.get_fpr();
		bf_occ = hits_filter.get_occupancy();
	} else if (args.seeds.empty()) {
		btllib::KmerBloomFilter bf(dbf_size, args.num_hashes, args.kmer_length);
		nthits::HitTable hits_table(hit_size);
		POPULATE_KMERS(hits_table)
		hits_table.save(args.out_file, args.thresh_min);
	} else {
		btllib::SeedBloomFilter bf(dbf_size, args.kmer_length, args.seeds, args.num_hashes);
		nthits::HitTable hits_table(hit_size);
		POPULATE_SEEDS(hits_table)
		hits_table.save(args.out_file, args.thresh_min);
	}
	timer.stop();
	timer.print_done();
	if (args.verbosity > 0) {
		print_bloom_filter_stats(bf_fpr, args.fpr, bf_occ);
	}

	return 0;
}
