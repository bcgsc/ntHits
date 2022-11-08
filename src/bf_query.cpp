#define PROGRAM_NAME "ntHits-BFQ"
#define PROGRAM_VERSION "0.0.1"
#define PROGRAM_DESCRIPTION "Query tool for ntHits' output Bloom filter"
#define PROGRAM_COPYRIGHT "Copyright 2022 Canada's Michael Smith Genome Science Centre"

#include <argparse/argparse.hpp>
#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <iostream>
#include <vector>

#include "nthits.hpp"
#include "user_interface.hpp"

struct ProgramArguments
{
	std::string bf_path;
	bool is_cbf;
	bool silent;
	std::vector<std::string> seeds;
};

ProgramArguments
parse_args(int argc, char** argv)
{
	ProgramArguments args;

	argparse::ArgumentParser parser(PROGRAM_NAME, PROGRAM_VERSION);
	parser.add_description(PROGRAM_DESCRIPTION);
	parser.add_epilog(PROGRAM_COPYRIGHT);

	parser.add_argument("bf_path").help("Input Bloom filter file").required();

	parser.add_argument("--cbf")
	    .help("Treat input file as a counting Bloom filter and output k-mer counts")
	    .default_value(false)
	    .implicit_value(true);

	parser.add_argument("-s", "--seeds")
	    .help("Spaced seed patterns separated with commas (e.g. 10101,11011)");

	parser.add_argument("--silent")
	    .help("Don't print logs to stdout")
	    .default_value(false)
	    .implicit_value(true);

	try {
		parser.parse_args(argc, argv);
	} catch (const std::runtime_error& err) {
		std::cerr << err.what() << std::endl;
		std::cerr << parser;
		std::exit(1);
	}

	args.bf_path = parser.get("bf_path");
	args.is_cbf = parser.get<bool>("--cbf");
	args.silent = parser.get<bool>("--silent");

	if (parser.is_used("-s")) {
		std::istringstream ss(parser.get("-s"));
		std::string seed;
		while (std::getline(ss, seed, ',')) {
			args.seeds.push_back(seed);
		}
	}

	return args;
}

#define GET_HASHES                                                                                 \
	const uint64_t* hashes;                                                                        \
	if (args.seeds.size() == 0) {                                                                  \
		btllib::NtHash nthash(kmer, bf.get_hash_num(), kmer.size());                               \
		nthash.roll();                                                                             \
		hashes = nthash.hashes();                                                                  \
	} else {                                                                                       \
		btllib::SeedNtHash nthash(kmer, args.seeds, bf.get_hash_num(), kmer.size());               \
		nthash.roll();                                                                             \
		hashes = nthash.hashes();                                                                  \
	}

int
main(int argc, char** argv)
{
	auto args = parse_args(argc, argv);

    omp_set_num_threads(1);

	if (!args.silent) {
		print_logo();
	}

	Timer timer;
	std::string kmer;

	if (!args.silent && args.is_cbf) {
		TIME_EXECUTION("Loading Bloom filter",
		               timer,
		               btllib::CountingBloomFilter<nthits::cbf_counter_t> bf(args.bf_path);)
		print_bloom_filter_stats(bf.get_fpr(), 1.0, bf.get_occupancy());
		std::cout << std::endl << "> " << std::flush;
		while (std::cin >> kmer) {
			GET_HASHES
			std::cout << (unsigned)bf.contains(hashes) << std::endl;
			std::cout << "> " << std::flush;
		}
	}
	if (args.silent && args.is_cbf) {
		btllib::CountingBloomFilter<nthits::cbf_counter_t> bf(args.bf_path);
		while (std::cin >> kmer) {
			GET_HASHES
			std::cout << (unsigned)bf.contains(hashes) << std::endl;
		}
	}
	if (!args.silent && !args.is_cbf) {
		TIME_EXECUTION("Loading Bloom filter", timer, btllib::BloomFilter bf(args.bf_path);)
		print_bloom_filter_stats(bf.get_fpr(), 1.0, bf.get_occupancy());
		std::cout << std::endl << "> " << std::flush;
		while (std::cin >> kmer) {
			GET_HASHES
			std::cout << (bf.contains(hashes) ? "TRUE" : "FALSE") << std::endl;
			std::cout << "> " << std::flush;
		}
	}
	if (args.silent && !args.is_cbf) {
		btllib::BloomFilter bf(args.bf_path);
		while (std::cin >> kmer) {
			GET_HASHES
			std::cout << (bf.contains(hashes) ? "TRUE" : "FALSE") << std::endl;
		}
	}

	return 0;
}