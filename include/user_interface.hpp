#ifndef NTHITS_UI_HPP
#define NTHITS_UI_HPP

#include <chrono>
#include <iostream>

#include "args.hpp"

#define LOGO                                                                                       \
	"       _          _ _              \n"                                                        \
	" _ __ | |_  /\\  /(_) |_ ___       \n"                                                        \
	"| '_ \\| __|/ /_/ / | __/ __|      \n"                                                        \
	"| | | | |_/ __  /| | |_\\__ \\     \n"                                                        \
	"|_| |_|\\__\\/ /_/ |_|\\__|___/    "

class Timer
{
  private:
	std::chrono::time_point<std::chrono::system_clock> t_start;
	std::chrono::time_point<std::chrono::system_clock> t_end;

  public:
	/**
	 * Register the current time as the timer's starting point.
	 */
	void start() { this->t_start = std::chrono::system_clock::now(); }

	/**
	 * Register the current time as the timer's finish point.
	 */
	void stop() { this->t_end = std::chrono::system_clock::now(); }

	/**
	 * Compute the difference between the start and stop points in seconds.
	 */
	[[nodiscard]] long double elapsed_seconds() const
	{
		std::chrono::duration<double> elapsed = t_end - t_start;
		return elapsed.count();
	}

	/**
	 * Get a human-readable representation of the elapsed time.
	 */
	[[nodiscard]] std::string to_string() const
	{
		return std::to_string(this->elapsed_seconds()) + "s";
	}

	/**
	 * Print a line-ender for logging.
	 */
	void print_done() const
	{
		std::cout << "\033[1;32m"
		          << "DONE"
		          << "\033[0m"
		          << "(" << this->to_string() << ")" << std::endl;
	}
};

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
	std::cout << "Sequence reading mode    : " << (args.long_mode ? "LONG" : "SHORT") << std::endl;
	if (args.seeds.size() > 0) {
		std::cout << "[-s] Spaced seed patterns:" << std::endl;
		for (const auto& seed : args.seeds) {
			std::cout << "  - " << seed << std::endl;
		}
		std::cout << "[-h] Hashes per seed     : " << args.num_hashes << std::endl;
	} else {
		std::cout << "[-k] k-mer length        : " << args.kmer_length << std::endl;
		std::cout << "[-h] Hashes per k-mer    : " << args.num_hashes << std::endl;
	}
	if (args.thresh_min > 0) {
		std::cout << "[-c] Filter threshold    : " << args.thresh_min << std::endl;
	}
	if (args.fpr > 0 && args.out_bloom) {
		std::cout << "[-p] Bloom filter FPR    : " << args.fpr << std::endl;
	}
	std::cout << "[-t] Number of threads   : " << args.num_threads << std::endl;
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
		std::cout << "- Output Bloom filter size  : " << hit_size << " bytes" << std::endl;
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

#endif // NTHITS_UI_HPP