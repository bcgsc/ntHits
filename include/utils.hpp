#ifndef NTHITS_UTILS_HPP
#define NTHITS_UTILS_HPP

#include <stdint.h>
#include <string>
#include <vector>

namespace nthits {

static const unsigned char b2r[256] = { 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 0
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 1
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 2
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 3
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', // 4   'A' 'C' 'G'
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', // 5   'T'
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', // 6   'a' 'c' 'g'
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', // 7   't'
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 8
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 9
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 10
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 11
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 12
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 13
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 14
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 15
	                                    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N' };

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

void
get_thresholds(std::vector<uint64_t> histogram, bool solid, size_t& hit_count, unsigned& hit_cap)
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

	if (hit_cap == 0 && solid) {
		hit_cap = min_cov;
	} else if (hit_cap == 0) {
		hit_cap = 1.75 * max_cov;
	}

	hit_count = histogram[1];
	for (unsigned i = 2; i <= hit_cap + 1; i++)
		hit_count -= histogram[i];
}

size_t
calc_bf_size(double num_elements, double num_hashes, int num_seeds, double fpr)
{
	double r = -num_hashes / log(1.0 - exp(log(fpr) / num_hashes));
	size_t m = ceil(num_elements * r);
	return m * std::max(1, num_seeds);
}

}

#endif // NTHITS_UTILS_HPP
