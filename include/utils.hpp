#ifndef NTHITS_UTILS_HPP
#define NTHITS_UTILS_HPP

#include <string>

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
getCanon(std::string& bMer);

#endif // NTHITS_UTILS_HPP
