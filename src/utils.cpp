#include "utils.hpp"

void
getCanon(std::string& bMer)
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