#include "utils.hpp"

void
nthits::to_canonical(std::string& bMer)
{
	int p = 0, hLen = (bMer.length() - 1) / 2;
	while (bMer[p] == b2r[(unsigned char)bMer[bMer.length() - 1 - p]]) {
		++p;
		if (p >= hLen)
			break;
	}
	if (bMer[p] > b2r[(unsigned char)bMer[bMer.length() - 1 - p]]) {
		for (int lIndex = p, rIndex = bMer.length() - 1 - p; lIndex <= rIndex; ++lIndex, --rIndex) {
			char tmp = b2r[(unsigned char)bMer[rIndex]];
			bMer[rIndex] = b2r[(unsigned char)bMer[lIndex]];
			bMer[lIndex] = tmp;
		}
	}
}

void
nthits::get_thresholds(
    std::vector<uint64_t> histogram,
    bool solid,
    size_t& hit_count,
    unsigned& hit_cap)
{
	int hist_index = 2, err_cov = 1;
	while (hist_index <= 10000 && histogram[hist_index] > histogram[hist_index + 1])
		hist_index++;

	err_cov = hist_index > 300 ? 1 : hist_index - 1;

	unsigned max_cov = err_cov;
	for (unsigned i = err_cov; i < 10002; i++) {
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
