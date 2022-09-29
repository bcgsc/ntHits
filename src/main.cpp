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

#define POPULATE(HITS_CONTAINER)                                                                   \
	for (const auto& file_path : args.input_files) {                                               \
		btllib::SeqReader reader(file_path, seq_reader_mode);                                      \
		_Pragma("omp parallel shared(reader)");                                                    \
		for (const auto& record : reader) {                                                        \
			nthits::process(record.seq, args.hit_cap, bf, cbf, HITS_CONTAINER);                    \
		}                                                                                          \
	}

int
main(int argc, char** argv)
{
	double sTime = omp_get_wtime();

	auto args = parse_arguments(argc, argv);
	unsigned dbf_size, cbf_size, hit_size;
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

	if (args.use_ntcard) {
		ntcard::NtCard ntc(args.kmer_length);
		for (const auto& file : args.input_files) {
			btllib::SeqReader reader(file, seq_reader_mode);
			for (const auto& record : reader) {
				ntc.process(record.seq);
			}
		}
		auto histArray = ntc.get_histogram(10000);

		int histIndex = 2, errCov = 1;
		while (histIndex <= 10000 && histArray[histIndex] > histArray[histIndex + 1])
			histIndex++;

		errCov = histIndex > 300 ? 1 : histIndex - 1;

		unsigned maxCov = errCov;
		for (unsigned i = errCov; i < 10002; i++) {
			if (histArray[i] >= histArray[maxCov])
				maxCov = i;
		}
		maxCov--;

		unsigned minCov = errCov;
		for (unsigned i = errCov; i < maxCov; i++) {
			if (histArray[i] <= histArray[minCov])
				minCov = i;
		}
		minCov--;

		if (args.solid) {
			if (args.hit_cap == 0)
				args.hit_cap = minCov;
			std::cerr << "Errors k-mer coverage: " << args.hit_cap << std::endl;
		} else {
			if (args.hit_cap == 0)
				args.hit_cap = 1.75 * maxCov;
			std::cerr << "Errors k-mer coverage: " << minCov << std::endl;
			std::cerr << "Median k-mer coverage: " << maxCov << std::endl;
			std::cerr << "Repeat k-mer coverage: " << args.hit_cap << std::endl;
		}

		dbf_size = args.bits * histArray[1];
		cbf_size = args.bytes * (histArray[1] - histArray[2]);
		size_t hitCount = histArray[1];
		for (unsigned i = 2; i <= args.hit_cap + 1; i++)
			hitCount -= histArray[i];
		hit_size = args.out_bloom ? hitCount * args.m : hitCount * 3;

		std::cerr << "Approximate# of distinct k-mers: " << histArray[1] << "\n";
		std::cerr << "Approximate# of solid k-mers: " << hitCount << "\n";
	} else {
		dbf_size = args.bits * args.f0;
		cbf_size = args.bytes * (args.f0 - args.f1);
		hit_size = args.m * args.fr;
		std::cerr << "Approximate# of distinct k-mers: " << args.f0 << "\n";
		std::cerr << "Approximate# of solid k-mers: " << args.fr << "\n";
	}

	btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);

	if (args.out_bloom && args.seeds.empty()) {
		btllib::KmerBloomFilter bf(dbf_size / 8, 3, args.kmer_length);
		btllib::KmerBloomFilter hits_filter(hit_size / 8, args.num_hashes + 1, args.kmer_length);
		POPULATE(hits_filter)
		hits_filter.save(out_file_prefix + ".bf");
	} else if (args.out_bloom) {
		btllib::SeedBloomFilter bf(dbf_size / 8, args.kmer_length, args.seeds, 3);
		btllib::SeedBloomFilter hits_filter(
		    hit_size / 8, args.kmer_length, args.seeds, args.num_hashes + 1);
		POPULATE(hits_filter)
		hits_filter.save(out_file_prefix + ".bf");
	} else if (args.seeds.empty()) {
		btllib::KmerBloomFilter bf(dbf_size / 8, 3, args.kmer_length);
		nthits::HitTable hits_table(hit_size);
		POPULATE(hits_table)
		hits_table.save(out_file_prefix + ".rep", args.hit_cap);
	} else {
		btllib::SeedBloomFilter bf(dbf_size / 8, args.kmer_length, args.seeds, 3);
		nthits::HitTable hits_table(hit_size);
		POPULATE(hits_table)
		hits_table.save(out_file_prefix + ".rep", args.hit_cap);
	}

	std::cerr << "Total time for computing repeat content in (sec): " << std::setprecision(4)
	          << std::fixed << omp_get_wtime() - sTime << "\n";
	return 0;
}
