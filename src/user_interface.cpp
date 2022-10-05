#include "user_interface.hpp"

#include <iostream>

void
Timer::start()
{
	this->t_start = clock();
}

void
Timer::stop()
{
	this->t_end = clock();
}

long double
Timer::elapsed_seconds() const
{
	return (long double)(this->t_end - this->t_start) / CLOCKS_PER_SEC;
}

std::string
Timer::to_string() const
{
	return std::to_string(this->elapsed_seconds()) + "s";
}

void
Timer::print_done() const
{
	std::cout << "\033[1;32m"
	          << "DONE"
	          << "\033[0m"
	          << "(" << this->to_string() << ")" << std::endl;
}

void
print_logo()
{
	std::cout << LOGO << "\tv" << PROGRAM_VERSION << std::endl << std::endl;
}

void
print_args(const ProgramArguments& args)
{
	std::cout << "Input files:" << std::endl;
	for (const auto& file : args.input_files) {
		std::cout << "  - " << file << std::endl;
	}
	if (args.seeds.size() > 0) {
		std::cout << "Spaced seed patterns (-s):" << std::endl;
		for (const auto& seed : args.seeds) {
			std::cout << "  - " << seed << std::endl;
		}
		std::cout << "Hashes per seed          (-h): " << args.num_hashes << std::endl;
	} else {
		std::cout << "k-mer length             (-k): " << args.kmer_length << std::endl;
		std::cout << "Hashes per k-mer         (-h): " << args.num_hashes << std::endl;
	}
	if (args.hit_cap > 0) {
		std::cout << "Filter threshold         (-c): " << args.hit_cap << std::endl;
	}
	std::cout << "Number of threads        (-t): " << args.num_threads << std::endl;
	std::cout << "Optimize file reading for long sequences (--long-mode): "
	          << (args.long_mode ? "YES" : "NO") << std::endl;
	std::cout << "Output mode ";
	if (args.out_bloom && args.solid) {
		std::cout << "(\033[1;34m--out-bloom\033[0m/\033[1;34m--solid\033[0m): ";
		std::cout << "Non-erroneous k-mers in a Bloom filter" << std::endl;
	} else if (args.out_bloom) {
		std::cout << "(\033[1;34m--out-bloom\033[0m/--solid): ";
		std::cout << "Repeated k-mers in a Bloom filter" << std::endl;
	} else if (args.solid) {
		std::cout << "(--out-bloom/\033[1;34m--solid\033[0m): ";
		std::cout << "Non-erroneous k-mers and counts in a table" << std::endl;
	} else {
		std::cout << "(--out-bloom/--solid): ";
		std::cout << "Repeated k-mers and counts in a table" << std::endl;
	}
	std::cout << std::endl;
}

void
print_ntcard_results(
    const size_t hit_count,
    const uint64_t num_distinct,
    const unsigned hit_cap,
    const bool hit_cap_changed)
{
	std::cout << "  - Estimated number of distinct k-mers  : " << num_distinct << std::endl;
	std::cout << "  - Estimated number of solid k-mers     : " << hit_count << std::endl;
	if (hit_cap_changed) {
		std::cout << "  - Estimated error k-mer threshold (-c) : " << hit_cap << std::endl;
	}
	std::cout << std::endl;
}

void
print_bloom_filter_stats(const double fpr, const double target_fpr, const double occupancy)
{
	if (fpr > target_fpr * 25) {
		std::cout << "  - Actual false positive rate (FPR): \033[1;31m" << fpr << "\033[0m"
		          << std::endl;
	} else if (fpr > target_fpr * 10) {
		std::cout << "  - Actual false positive rate (FPR): \033[1;33m" << fpr << "\033[0m"
		          << std::endl;
	} else {
		std::cout << "  - Actual false positive rate (FPR): \033[1;32m" << fpr << "\033[0m"
		          << std::endl;
	}
	std::cout << "  - Bloom filter occupancy: " << occupancy << std::endl;
}