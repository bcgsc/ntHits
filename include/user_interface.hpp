#ifndef NTHITS_UI_HPP
#define NTHITS_UI_HPP

#include <chrono>
#include <iostream>
#include <thread>

#include "human_readable_strings.hpp"
#include "logo.hpp"
#include "nthits.hpp"
#include "version.hpp"

enum Color
{
  FG_RED = 31,
  FG_GREEN = 32,
  FG_BLUE = 34,
  FG_YELLOW = 33,
  FG_DEFAULT = 39,
  BG_RED = 41,
  BG_GREEN = 42,
  BG_BLUE = 44,
  BG_YELLOW = 43,
  BG_DEFAULT = 49
};

inline std::string
color(const std::string& text, Color color)
{
  return "\033[1;" + std::to_string(color) + "m" + text + "\033[0m";
}

inline std::string
comma_sep(uint64_t val)
{
  std::string val_str = std::to_string(val);
  unsigned i = val_str.size() % 3;
  std::string result = i > 0 ? val_str.substr(0, i) + "," : "";
  for (; i + 3 <= val_str.size(); i += 3) {
    result += val_str.substr(i, 3) + ",";
  }
  return result.substr(0, result.size() - 1);
}

class Timer
{
private:
  std::chrono::time_point<std::chrono::system_clock> t_start;

public:
  void start(std::string message)
  {
    std::cout << message << "... " << std::flush;
    t_start = std::chrono::system_clock::now();
  }

  void stop()
  {
    auto t_end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = t_end - t_start;
    std::cout << color("DONE", Color::FG_GREEN) << " (" << elapsed.count() << "s)" << std::endl;
  }
};

void
print_logo()
{
  std::cerr << LOGO << "\tv" << nthits::VERSION << std::endl << std::endl;
}

void
print_updated_params(uint64_t kmer_count,
                     size_t hit_count,
                     uint64_t num_distinct,
                     size_t bf_size,
                     unsigned hit_cap,
                     bool hit_cap_changed,
                     size_t cbf_size,
                     bool out_bloom,
                     size_t hit_size,
                     size_t ex_size,
                     unsigned verbosity)
{
  std::cout << "- Total number of k-mers                  : " << comma_sep(kmer_count) << std::endl;
  std::cout << "- Number of distinct k-mers               : " << comma_sep(num_distinct)
            << std::endl;
  std::cout << "- Number of filtered k-mers               : " << comma_sep(hit_count) << std::endl;
  if (hit_cap_changed) {
    std::cout << "- Estimated error k-mer threshold         : " << hit_cap << std::endl;
  }
  if (verbosity > 1) {
    if (bf_size > 1)
      std::cout << "- Distinct k-mers Bloom filter size       : " << comma_sep(bf_size)
                << std::endl;
    if (cbf_size > 1)
      std::cout << "- Counting Bloom filter size              : " << comma_sep(cbf_size)
                << std::endl;
  }
  if (out_bloom) {
    std::cout << "- Output Bloom filter size                : " << comma_sep(hit_size) << std::endl;
    if (ex_size > 0) {
      std::cout << "- Excluded k-mers Bloom filter size       : " << comma_sep(ex_size)
                << std::endl;
    }
  }
  std::cout << std::endl;
}

void
print_bloom_filter_stats(const double fpr, const double target_fpr, const double occupancy)
{
  std::string fpr_str;
  if (fpr > target_fpr * 25) {
    fpr_str = color(std::to_string(fpr), Color::FG_RED);
  } else if (fpr > target_fpr * 10) {
    fpr_str = color(std::to_string(fpr), Color::FG_YELLOW);
  } else {
    fpr_str = color(std::to_string(fpr), Color::FG_GREEN);
  }
  std::cout << "  - Actual false positive rate (FPR) : " << fpr_str << std::endl;
  std::cout << "  - Bloom filter occupancy           : " << occupancy << std::endl;
}

#endif // NTHITS_UI_HPP