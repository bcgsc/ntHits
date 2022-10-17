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
	double fpr;
	unsigned thresh_min = 0;
	std::string out_file;
	unsigned verbosity;
	bool out_bloom;
	bool solid;
	bool long_mode;
	std::string histogram_path;
	std::vector<std::string> input_files;
	std::vector<std::string> seeds;
};

ProgramArguments
parse_arguments(int argc, char** argv);

#endif // NTHITS_ARGS_HPP
