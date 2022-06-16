#include "args.hpp"

#include <argparse/argparse.hpp>
#include <sstream>

void
ProgramArguments::parse(int argc, char** argv)
{
	auto default_args = argparse::default_arguments::none;
	argparse::ArgumentParser parser(PROGRAM_NAME, PROGRAM_VERSION, default_args);
	parser.add_description(PROGRAM_DESCRIPTION);
	parser.add_epilog(PROGRAM_COPYRIGHT);

	parser.add_argument("-t", "--threads")
	    .help("Number of parallel threads")
	    .default_value(16U)
	    .scan<'u', unsigned>();

	parser.add_argument("-k", "--kmer")
	    .help("k-mer length")
	    .default_value(64U)
	    .scan<'u', unsigned>();

	parser.add_argument("-h", "--hashes")
	    .help("Number of hashes to generate per k-mer/spaced seed")
	    .default_value(4U)
	    .scan<'u', unsigned>();

	parser.add_argument("-c", "--cutoff")
	    .help("k-mer cutoff threshold")
	    .required()
	    .scan<'u', unsigned>();

	parser.add_argument("-p", "--prefix")
	    .help("Output files' prefix")
	    .default_value(std::string("repeat"));

	parser.add_argument("-s", "--seeds")
	    .help("If specified, use given spaced seeds separated by commas (e.g. 10101,11011)");

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

	t = parser.get<unsigned>("-t");
	k = parser.get<unsigned>("-k");
	h = parser.get<unsigned>("-h");
	m = parser.get<unsigned>("-b");
	hit_cap = parser.get<unsigned>("-c");
	prefix = parser.get("-p");
	_out_bloom = parser.get<bool>("--outbloom");
	_solid = parser.get<bool>("--solid");
	_long_mode = parser.get<bool>("--long-mode");

	_use_ntcard = true;

	if (parser.is_used("-F")) {
		F0 = parser.get<unsigned>("-F");
		_use_ntcard = false;
	}

	if (parser.is_used("-r")) {
		fr = parser.get<unsigned>("-r");
		_use_ntcard = false;
	}

	if (parser.is_used("-f")) {
		f1 = parser.get<unsigned>("-f");
		_use_ntcard = false;
	}

	if (parser.is_used("-s")) {
		std::istringstream ss(parser.get("-s"));
		std::string seed;
		while (std::getline(ss, seed, ',')) {
			seeds.push_back(seed);
		}
	}

	try {
		input_files = parser.get<std::vector<std::string> >("files");
	} catch (std::logic_error& e) {
		std::cerr << "No files provided" << std::endl;
	}
}