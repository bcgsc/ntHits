#include "args.hpp"

#include <argparse/argparse.hpp>
#include <sstream>

ProgramArguments
parse_arguments(int argc, char** argv)
{
	ProgramArguments args;

	auto default_args = argparse::default_arguments::none;
	argparse::ArgumentParser parser(PROGRAM_NAME, PROGRAM_VERSION, default_args);
	parser.add_description(PROGRAM_DESCRIPTION);
	parser.add_epilog(PROGRAM_COPYRIGHT);

	parser.add_argument("-t", "--threads")
	    .help("Number of parallel threads")
	    .default_value(16U)
	    .scan<'u', unsigned>();

	parser.add_argument("-k", "--kmer")
	    .help("k-mer length, ignored if using spaced seeds (-s)")
	    .default_value(64U)
	    .scan<'u', unsigned>();

	parser.add_argument("-h", "--hashes")
	    .help("Number of hashes to generate per k-mer/spaced seed")
	    .default_value(4U)
	    .scan<'u', unsigned>();

	parser.add_argument("-c", "--cutoff")
	    .help("k-mer cutoff threshold")
	    .default_value(0U)
	    .scan<'u', unsigned>();

	parser.add_argument("-p", "--prefix")
	    .help("Output files' prefix")
	    .default_value(std::string("repeat"));

	parser.add_argument("-s", "--seeds")
	    .help("If specified, use spaced seeds (separate with commas, e.g. 10101,11011)");

	parser.add_argument("--outbloom")
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

	parser.add_argument("-b", "--bit").default_value(16U).scan<'u', unsigned>();
	parser.add_argument("-F").scan<'u', unsigned>();
	parser.add_argument("-f").scan<'u', unsigned>();
	parser.add_argument("-r").scan<'u', unsigned>();

	parser.add_argument("files").help("Input files").required().remaining();

	try {
		parser.parse_args(argc, argv);
	} catch (const std::runtime_error& err) {
		std::cerr << err.what() << std::endl;
		std::cerr << parser;
		std::exit(1);
	}

	args.num_threads = parser.get<unsigned>("-t");
	args.kmer_length = parser.get<unsigned>("-k");
	args.num_hashes = parser.get<unsigned>("-h");
	args.m = parser.get<unsigned>("-b");
	args.hit_cap = parser.get<unsigned>("-c");
	args.prefix = parser.get("-p");
	args.out_bloom = parser.get<bool>("--outbloom");
	args.solid = parser.get<bool>("--solid");
	args.long_mode = parser.get<bool>("--long-mode");

	args.use_ntcard = true;

	if (parser.is_used("-F")) {
		args.f0 = parser.get<unsigned>("-F");
		args.use_ntcard = false;
	}

	if (parser.is_used("-r")) {
		args.fr = parser.get<unsigned>("-r");
		args.use_ntcard = false;
	}

	if (parser.is_used("-f")) {
		args.f1 = parser.get<unsigned>("-f");
		args.use_ntcard = false;
	}

	if (parser.is_used("-s")) {
		std::istringstream ss(parser.get("-s"));
		std::string seed;
		while (std::getline(ss, seed, ',')) {
			args.seeds.push_back(seed);
			args.kmer_length = seed.size();
		}
	}

	try {
		args.input_files = parser.get<std::vector<std::string> >("files");
	} catch (std::logic_error& e) {
		std::cerr << "No files provided" << std::endl;
	}

	return args;
}