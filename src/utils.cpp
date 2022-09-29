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

inline void
nthits::get_thresholds(
    std::vector<uint64_t> histogram,
    bool solid,
    size_t& hit_count,
    unsigned& hit_cap)
{
	unsigned min_cov, max_cov;
	int histIndex = 2, errCov = 1;
	while (histIndex <= 10000 && histogram[histIndex] > histogram[histIndex + 1])
		histIndex++;

	errCov = histIndex > 300 ? 1 : histIndex - 1;

	unsigned max_coverage = errCov;
	for (unsigned i = errCov; i < 10002; i++) {
		if (histogram[i] >= histogram[max_coverage])
			max_coverage = i;
	}
	max_coverage--;

	unsigned min_coverage = errCov;
	for (unsigned i = errCov; i < max_coverage; i++) {
		if (histogram[i] <= histogram[min_coverage])
			min_coverage = i;
	}
	min_coverage--;

	if (hit_cap == 0 && solid) {
		hit_cap = min_cov;
	} else if (hit_cap == 0) {
		hit_cap = 1.75 * max_cov;
	}

	hit_count = histogram[1];
	for (unsigned i = 2; i <= hit_cap + 1; i++)
		hit_count -= histogram[i];
}
