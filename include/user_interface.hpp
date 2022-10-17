#ifndef NTHITS_UI_HPP
#define NTHITS_UI_HPP

#include <chrono>

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
	std::clock_t t_start;
	std::clock_t t_end;

  public:
	/**
	 * Register the current time as the timer's starting point.
	 */
	void start();

	/**
	 * Register the current time as the timer's finish point.
	 */
	void stop();

	/**
	 * Compute the difference between the start and stop points in seconds.
	 */
	[[nodiscard]] long double elapsed_seconds() const;

	/**
	 * Get a human-readable representation of the elapsed time.
	 */
	[[nodiscard]] std::string to_string() const;

	/**
	 * Print a line-ender for logging.
	 */
	void print_done() const;
};

void
print_logo();

void
print_args(const ProgramArguments& args);

void
print_updated_params(
    size_t hit_count,
    uint64_t num_distinct,
    unsigned hit_cap,
    bool hit_cap_changed,
    bool out_bloom,
    size_t hit_size);

void
print_bloom_filter_stats(const double fpr, const double target_fpr, const double occupancy);

#endif // NTHITS_UI_HPP