#ifndef NTHITS_ARGS_HPP
#define NTHITS_ARGS_HPP

#define PROGRAM_NAME "ntHits"
#define PROGRAM_VERSION "0.0.1"
#define PROGRAM_DESCRIPTION "Filters k-mers based on counts (cmin <= count <= cmax) in input files"
#define PROGRAM_COPYRIGHT "Copyright 2022 Canada's Michael Smith Genome Science Centre"

#include <argparse/argparse.hpp>
#include <limits>
#include <string>
#include <vector>

#include "nthits.hpp"

const nthits::cbf_counter_t CBF_COUNTER_MAX = std::numeric_limits<nthits::cbf_counter_t>::max() - 1;

struct ProgramArguments
{
	unsigned num_threads;
	unsigned kmer_length;
	unsigned num_hashes;
	double fpr;
	unsigned min_count = 0, max_count = CBF_COUNTER_MAX;
	std::string out_file;
	unsigned verbosity;
	bool out_bloom;
	bool solid;
	bool long_mode;
	std::string histogram_path;
	std::vector<std::string> input_files;
	std::vector<std::string> seeds;

	ProgramArguments(int argc, char** argv);

	void print();
};

ProgramArguments::ProgramArguments(int argc, char** argv)
{
	auto default_args = argparse::default_arguments::none;
	argparse::ArgumentParser parser(PROGRAM_NAME, PROGRAM_VERSION, default_args);
	parser.add_description(PROGRAM_DESCRIPTION);
	parser.add_epilog(PROGRAM_COPYRIGHT);

	parser.add_argument("-t", "--threads")
	    .help("Number of parallel threads")
	    .default_value(4U)
	    .scan<'u', unsigned>();

	parser.add_argument("-k", "--kmer")
	    .help("k-mer length, ignored if using spaced seeds (-s)")
	    .default_value(64U)
	    .scan<'u', unsigned>();

	parser.add_argument("-h", "--hashes")
	    .help("Number of hashes to generate per k-mer/spaced seed")
	    .default_value(3U)
	    .scan<'u', unsigned>();

	parser.add_argument("-f", "--frequencies")
	    .help("Frequency histogram file (e.g. from ntCard)")
	    .required();

	parser.add_argument("-p")
	    .help("Target Bloom filter false positive rate")
	    .default_value((double)0.0001)
	    .scan<'g', double>();

	parser.add_argument("-cmin", "--min-count")
	    .help("Minimum k-mer count (>1)")
	    .scan<'u', unsigned>();

	parser.add_argument("-cmax", "--max-count")
	    .help("Maximum k-mer count (<" + std::to_string(max_count + 1) + ")")
	    .scan<'u', unsigned>();

	parser.add_argument("-o", "--out").help("Output file's name").required();

	parser.add_argument("-s", "--seeds")
	    .help("If specified, use spaced seeds (separate with commas, e.g. 10101,11011)");

	parser.add_argument("--out-bloom")
	    .help("Output the most frequent k-mers in a Bloom filter")
	    .default_value(false)
	    .implicit_value(true);

	parser.add_argument("--solid")
	    .help("Output the solid k-mers (non-erroneous k-mers)")
	    .default_value(false)
	    .implicit_value(true);

	parser.add_argument("--long-mode")
	    .help("Optimize data reader for long sequences (>5kbp)")
	    .default_value(false)
	    .implicit_value(true);

	parser.add_argument("-v")
	    .action([&](const auto&) { ++verbosity; })
	    .append()
	    .help("Level of details printed to stdout (-v: normal, -vv detailed)")
	    .default_value(false)
	    .implicit_value(true)
	    .nargs(0);

	parser.add_argument("files").help("Input files").required().remaining();

	try {
		parser.parse_args(argc, argv);
	} catch (const std::runtime_error& err) {
		std::cerr << err.what() << std::endl;
		std::cerr << parser;
		std::exit(1);
	}

	if (parser.is_used("-cmin")) {
		min_count = parser.get<unsigned>("-cmin");
	}
	if (parser.is_used("-cmax")) {
		max_count = parser.get<unsigned>("-cmax");
	}

	if (min_count <= 1 || max_count > CBF_COUNTER_MAX) {
		std::cerr << "Invalid k-mer count range, must be cmin > 1 and cmax < " << CBF_COUNTER_MAX
		          << std::endl;
		std::exit(1);
	}

	num_threads = parser.get<unsigned>("-t");
	kmer_length = parser.get<unsigned>("-k");
	num_hashes = parser.get<unsigned>("-h");
	histogram_path = parser.get("-f");
	fpr = parser.get<double>("-p");
	out_file = parser.get("-o");
	out_bloom = parser.get<bool>("--out-bloom");
	solid = parser.get<bool>("--solid");
	long_mode = parser.get<bool>("--long-mode");

	if (parser.is_used("-s")) {
		std::istringstream ss(parser.get("-s"));
		std::string seed;
		while (std::getline(ss, seed, ',')) {
			seeds.push_back(seed);
			kmer_length = seed.size();
		}
	}

	try {
		input_files = parser.get<std::vector<std::string> >("files");
	} catch (std::logic_error& e) {
		std::cerr << "No files provided" << std::endl;
		std::cout << parser << std::endl;
		exit(0);
	}
}

void
ProgramArguments::print()
{
	std::cout << "Input files:" << std::endl;
	for (const auto& file : input_files) {
		std::cout << "  - " << file << std::endl;
	}
	std::cout << "Sequence reading mode       : " << (long_mode ? "LONG" : "SHORT") << std::endl;
	if (seeds.size() > 0) {
		std::cout << "[-s] Spaced seed patterns   :" << std::endl;
		for (const auto& seed : seeds) {
			std::cout << "  - " << seed << std::endl;
		}
		std::cout << "[-h] Hashes per seed        : " << num_hashes << std::endl;
	} else {
		std::cout << "[-k] k-mer length           : " << kmer_length << std::endl;
		std::cout << "[-h] Hashes per k-mer       : " << num_hashes << std::endl;
	}
	if (fpr > 0 && out_bloom) {
		std::cout << "[-p] Bloom filter FPR       : " << fpr << std::endl;
	}
	std::cout << "[-t] Number of threads      : " << num_threads << std::endl;
	if (min_count > 0) {
		std::cout << "[-cmin] Min. k-mer count    : " << min_count << std::endl;
	}
	if (max_count < CBF_COUNTER_MAX) {
		std::cout << "[-cmax] Max. k-mer count    : " << max_count << std::endl;
	}
	std::cout << std::endl;
}

#endif // NTHITS_ARGS_HPP
