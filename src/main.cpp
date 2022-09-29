#include "args.hpp"
#include "ntcard.hpp"
#include "nthits.hpp"
#include "user_interface.hpp"
#include "utils.hpp"

#include <argparse/argparse.hpp>
#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <btllib/nthash.hpp>
#include <btllib/seq_reader.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

#define POPULATE_KMERS(HITS_CONTAINER)                                                             \
	for (const auto& file_path : args.input_files) {                                               \
		btllib::SeqReader reader(file_path, seq_reader_mode);                                      \
		_Pragma("omp parallel shared(reader)");                                                    \
		for (const auto& record : reader) {                                                        \
			nthits::process(record.seq, args.hit_cap, bf, cbf, HITS_CONTAINER);                    \
		}                                                                                          \
	}

#define POPULATE_SEEDS(HITS_CONTAINER)                                                             \
	for (const auto& file_path : args.input_files) {                                               \
		btllib::SeqReader reader(file_path, seq_reader_mode);                                      \
		_Pragma("omp parallel shared(reader)");                                                    \
		for (const auto& record : reader) {                                                        \
			for (const auto& seed : args.seeds) {                                                  \
				nthits::process(record.seq, seed, args.hit_cap, bf, cbf, HITS_CONTAINER);          \
			}                                                                                      \
		}                                                                                          \
	}

int
main(int argc, char** argv)
{
	print_logo();

	auto args = parse_arguments(argc, argv);
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

	Timer timer;

	unsigned dbf_size, cbf_size, hit_size;
	if (args.use_ntcard) {
		std::cout << "Running ntCard... " << std::flush;
		timer.start();
		ntcard::NtCard ntc(args.kmer_length);
		for (const auto& file : args.input_files) {
			btllib::SeqReader reader(file, seq_reader_mode);
			for (const auto& record : reader) {
				ntc.process(record.seq);
			}
		}
		auto hist = ntc.get_histogram(10000);
		timer.stop();
		timer.print_done();

		size_t hit_count;
		unsigned given_hit_cap = args.hit_cap;
		nthits::get_thresholds(hist, args.solid, hit_count, args.hit_cap);
		bool hit_cap_changed = given_hit_cap != args.hit_cap;

		if (args.verbosity > 0) {
			print_ntcard_results(hit_count, hist[1], args.hit_cap, hit_cap_changed);
		}

		dbf_size = args.bits * hist[1];
		cbf_size = args.bytes * (hist[1] - hist[2]);
		hit_size = args.out_bloom ? hit_count * args.m : hit_count * 3;
	} else {
		dbf_size = args.bits * args.f0;
		cbf_size = args.bytes * (args.f0 - args.f1);
		hit_size = args.m * args.fr;
	}

	std::cout << "Running ntHits... " << std::flush;
	timer.start();
	btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
	double bf_fpr = 0, bf_occ = 0;
	if (args.out_bloom && args.seeds.empty()) {
		btllib::KmerBloomFilter bf(dbf_size / 8, args.num_hashes, args.kmer_length);
		btllib::KmerBloomFilter hits_filter(hit_size / 8, args.num_hashes, args.kmer_length);
		POPULATE_KMERS(hits_filter)
		hits_filter.save(args.out_file);
		bf_fpr = hits_filter.get_fpr();
		bf_occ = hits_filter.get_occupancy();
	} else if (args.out_bloom) {
		btllib::SeedBloomFilter bf(dbf_size / 8, args.kmer_length, args.seeds, args.num_hashes);
		btllib::SeedBloomFilter hits_filter(
		    hit_size / 8, args.kmer_length, args.seeds, args.num_hashes);
		POPULATE_SEEDS(hits_filter)
		hits_filter.save(args.out_file);
		bf_fpr = hits_filter.get_fpr();
		bf_occ = hits_filter.get_occupancy();
	} else if (args.seeds.empty()) {
		btllib::KmerBloomFilter bf(dbf_size / 8, args.num_hashes, args.kmer_length);
		nthits::HitTable hits_table(hit_size);
		POPULATE_KMERS(hits_table)
		hits_table.save(args.out_file, args.hit_cap);
	} else {
		btllib::SeedBloomFilter bf(dbf_size / 8, args.kmer_length, args.seeds, args.num_hashes);
		nthits::HitTable hits_table(hit_size);
		POPULATE_SEEDS(hits_table)
		hits_table.save(args.out_file, args.hit_cap);
	}
	timer.stop();
	timer.print_done();
	if (args.verbosity > 0) {
		print_bloom_filter_stats(bf_fpr, bf_occ);
	}

	return 0;
}
