#ifndef NTHITS_ARGS_HPP
#define NTHITS_ARGS_HPP

#define PROGRAM_NAME "ntHits"
#define PROGRAM_VERSION "0.0.1"
#define PROGRAM_DESCRIPTION "Reports the most frequent k-mers in input files."
#define PROGRAM_COPYRIGHT "Copyright 2019 Canada's Michael Smith Genome Science Centre"

#include <string>
#include <vector>

struct ProgramArguments
{
	unsigned num_threads;
	unsigned kmer_length;
	unsigned num_hashes;
	unsigned hit_cap = 0;
	unsigned bits = 7;
	unsigned bytes = 6;
	unsigned m = 16;
	std::string prefix;
	unsigned f0, f1, fr;
	bool out_bloom;
	bool solid;
	bool long_mode;
	bool use_ntcard;
	std::vector<std::string> input_files;
	std::vector<std::string> seeds;
};

ProgramArguments
parse_arguments(int argc, char** argv);

#endif // NTHITS_ARGS_HPP
