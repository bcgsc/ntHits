#include "args.hpp"
#include "ntcard.hpp"
#include "nthits.hpp"
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
	auto args = parse_arguments(argc, argv);

	unsigned seq_reader_mode;
	if (args.long_mode) {
		seq_reader_mode = btllib::SeqReader::Flag::LONG_MODE;
	} else {
		seq_reader_mode = btllib::SeqReader::Flag::SHORT_MODE;
	}

	std::string out_file_prefix;
	if (args.out_file.empty() && args.solid) {
		out_file_prefix = "solids";
	} else if (args.out_file.empty() && !args.solid) {
		out_file_prefix = "repeat";
	} else {
		out_file_prefix = args.out_file;
	}

	omp_set_num_threads(args.num_threads);

	unsigned dbf_size, cbf_size, hit_size;
	if (args.use_ntcard) {
		ntcard::NtCard ntc(args.kmer_length);
		for (const auto& file : args.input_files) {
			btllib::SeqReader reader(file, seq_reader_mode);
			for (const auto& record : reader) {
				ntc.process(record.seq);
			}
		}
		auto hist = ntc.get_histogram(10000);

		size_t hit_count;
		nthits::get_thresholds(hist, args.solid, hit_count, args.hit_cap);

		dbf_size = args.bits * hist[1];
		cbf_size = args.bytes * (hist[1] - hist[2]);
		hit_size = args.out_bloom ? hit_count * args.m : hit_count * 3;
	} else {
		dbf_size = args.bits * args.f0;
		cbf_size = args.bytes * (args.f0 - args.f1);
		hit_size = args.m * args.fr;
	}

	btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
	if (args.out_bloom && args.seeds.empty()) {
		btllib::KmerBloomFilter bf(dbf_size / 8, args.num_hashes, args.kmer_length);
		btllib::KmerBloomFilter hits_filter(hit_size / 8, args.num_hashes, args.kmer_length);
		POPULATE_KMERS(hits_filter)
		hits_filter.save(out_file_prefix + ".bf");
	} else if (args.out_bloom) {
		btllib::SeedBloomFilter bf(dbf_size / 8, args.kmer_length, args.seeds, args.num_hashes);
		btllib::SeedBloomFilter hits_filter(
		    hit_size / 8, args.kmer_length, args.seeds, args.num_hashes);
		POPULATE_SEEDS(hits_filter)
		hits_filter.save(out_file_prefix + ".bf");
	} else if (args.seeds.empty()) {
		btllib::KmerBloomFilter bf(dbf_size / 8, args.num_hashes, args.kmer_length);
		nthits::HitTable hits_table(hit_size);
		POPULATE_KMERS(hits_table)
		hits_table.save(out_file_prefix + ".rep", args.hit_cap);
	} else {
		btllib::SeedBloomFilter bf(dbf_size / 8, args.kmer_length, args.seeds, args.num_hashes);
		nthits::HitTable hits_table(hit_size);
		POPULATE_SEEDS(hits_table)
		hits_table.save(out_file_prefix + ".rep", args.hit_cap);
	}

	return 0;
}
