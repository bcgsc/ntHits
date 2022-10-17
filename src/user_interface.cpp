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
	std::cout << "Sequence reading mode: " << (args.long_mode ? "LONG" : "SHORT") << std::endl;
	if (args.seeds.size() > 0) {
		std::cout << "[-s] Spaced seed patterns:" << std::endl;
		for (const auto& seed : args.seeds) {
			std::cout << "  - " << seed << std::endl;
		}
		std::cout << "[-h]   Hashes per seed   : " << args.num_hashes << std::endl;
	} else {
		std::cout << "[-k]   k-mer length      : " << args.kmer_length << std::endl;
		std::cout << "[-h]   Hashes per k-mer  : " << args.num_hashes << std::endl;
	}
	if (args.thresh_min > 0) {
		std::cout << "[-c]   Filter threshold  : " << args.thresh_min << std::endl;
	}
	if (args.fpr > 0 && args.out_bloom) {
		std::cout << "[-fpr] Bloom filter FPR  : " << args.fpr << std::endl;
	}
	std::cout << "[-t]   Number of threads : " << args.num_threads << std::endl;
	if (args.out_bloom && args.solid) {
		std::cout << "[\033[1;34m--out-bloom\033[0m/\033[1;34m--solid\033[0m]    : ";
		std::cout << "Non-erroneous k-mers in a Bloom filter" << std::endl;
	} else if (args.out_bloom) {
		std::cout << "[\033[1;34m--out-bloom\033[0m/--solid]    : ";
		std::cout << "Repeated k-mers in a Bloom filter" << std::endl;
	} else if (args.solid) {
		std::cout << "[--out-bloom/\033[1;34m--solid\033[0m]    : ";
		std::cout << "Non-erroneous k-mers and counts in a table" << std::endl;
	} else {
		std::cout << "[--out-bloom/--solid]    : ";
		std::cout << "Repeated k-mers and counts in a table" << std::endl;
	}
	std::cout << std::endl;
}

void
print_updated_params(
    size_t hit_count,
    uint64_t num_distinct,
    unsigned hit_cap,
    bool hit_cap_changed,
    bool out_bloom,
    size_t hit_size)
{
	std::cout << "- Number of distinct k-mers : " << num_distinct << std::endl;
	std::cout << "- Number of filtered k-mers : " << hit_count << std::endl;
	if (hit_cap_changed) {
		std::cout << "- Estimated error k-mer threshold : " << hit_cap << std::endl;
	}
	if (out_bloom) {
		std::cout << "- Output Bloom filter size : " << hit_size << std::endl;
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