#ifndef NTHITS_UTILS_HPP
#define NTHITS_UTILS_HPP

#include <cmath>
#include <fstream>
#include <stdint.h>
#include <string>
#include <vector>

namespace nthits {

static const unsigned char b2r[256] = {
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'A', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

void
to_canonical(std::string& bmer)
{
  int p = 0, h_len = (bmer.length() - 1) / 2;
  while (bmer[p] == b2r[(unsigned char)bmer[bmer.length() - 1 - p]]) {
    ++p;
    if (p >= h_len)
      break;
  }
  if (bmer[p] > b2r[(unsigned char)bmer[bmer.length() - 1 - p]]) {
    for (int l_idx = p, r_idx = bmer.length() - 1 - p; l_idx <= r_idx; ++l_idx, --r_idx) {
      char tmp = b2r[(unsigned char)bmer[r_idx]];
      bmer[r_idx] = b2r[(unsigned char)bmer[l_idx]];
      bmer[l_idx] = tmp;
    }
  }
}

inline std::vector<uint64_t>
load_ntcard_histogram(const std::string& path)
{
  std::vector<uint64_t> hist;
  std::ifstream hist_file(path);
  std::string freq;
  uint64_t value;
  while (hist_file >> freq >> value) {
    hist.push_back(value);
  }
  return hist;
}

void
get_thresholds(std::vector<uint64_t> histogram,
               bool solid,
               size_t& hit_count,
               size_t& ex_count,
               unsigned& min_count,
               unsigned max_count)
{
  int hist_index = 2, err_cov = 1;
  while (hist_index <= (int)histogram.size() - 2 &&
         histogram[hist_index] > histogram[hist_index + 1])
    hist_index++;

  err_cov = hist_index > 300 ? 1 : hist_index - 1;

  unsigned max_cov = err_cov;
  for (unsigned i = err_cov; i < histogram.size(); i++) {
    if (histogram[i] >= histogram[max_cov])
      max_cov = i;
  }
  max_cov--;

  unsigned min_cov = err_cov;
  for (unsigned i = err_cov; i < max_cov; i++) {
    if (histogram[i] <= histogram[min_cov])
      min_cov = i;
  }
  min_cov--;

  if (min_count == 0 && solid) {
    min_count = min_cov;
  } else if (min_count == 0) {
    min_count = 1.75 * max_cov;
  }

  if (max_count < histogram.size()) {
    hit_count = 0;
    for (unsigned i = min_count; i <= max_count; i++)
      hit_count += histogram[i];
  } else {
    hit_count = histogram[1];
    for (unsigned i = 2; i <= min_count + 1; i++)
      hit_count -= histogram[i];
  }

  ex_count = histogram[1];
  for (unsigned i = 2; i <= std::min(max_count, (unsigned)histogram.size() - 1); i++)
    ex_count -= histogram[i];
}

size_t
get_bf_size(double num_elements, double num_hashes, int num_seeds, double fpr)
{
  double r = -num_hashes / log(1.0 - exp(log(fpr) / num_hashes));
  size_t m = ceil(num_elements * r);
  return m * std::max(1, num_seeds);
}

}

#endif // NTHITS_UTILS_HPP
