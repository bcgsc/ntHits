#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "vendor/btl_bloomfilter/BloomFilter.hpp"
#include "CBFilter.hpp"
#include "ntcard.hpp"
#include "vendor/ntHash/ntHashIterator.hpp"

#include "Uncompress.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "nthits"

static const char VERSION_MESSAGE[] =
    PROGRAM " 0.1.1 \n"
            "Written by Hamid Mohamadi.\n"
            "Copyright 2019 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... FILES...\n"
    "Reports the most frequent k-mers in FILES(>=1).\n"
    "Accepatble file formats: fastq, fasta, gz, bz, zip.\n"
    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N	use N parallel threads [16]\n"
    "  -k, --kmer=N	the length of kmer [64]\n"
    "  -c, --cutoff=N	the maximum coverage of kmer in output\n"
    "  -p, --pref=STRING	the prefix for output file name [repeat]\n"
    "  --outbloom	output the most frequent k-mers in a Bloom filter\n"
    "  --solid	output the solid k-mers (non-errornous k-mers)\n"
    "  --help	display this help and exit\n"
    "  --version	output version information and exit\n"
    "\n"
    "Report bugs to hmohamadi@bcgsc.ca.\n";

using namespace std;

namespace opt {
size_t t = 16;
size_t k = 64;
size_t h = 4;
size_t hitCap = 0;
size_t bytes = 6;
size_t bits = 7;
size_t m = 16;
size_t dbfSize;
size_t cbfSize;
size_t hitSize;
string prefix;
size_t F0;
size_t f1;
size_t fr;
int outbloom = 0;
int solid = 0;
int eval = 0;
} // namespace opt

static const char shortopts[] = "t:h:k:b:p:r:c:F:f:";

enum
{
	OPT_HELP = 1,
	OPT_VERSION
};

static const struct option longopts[] = { { "threads", required_argument, NULL, 't' },
	                                      { "kmer", required_argument, NULL, 'k' },
	                                      { "hash", required_argument, NULL, 'h' },
	                                      { "bit", required_argument, NULL, 'b' },
	                                      { "cutoff", required_argument, NULL, 'c' },
	                                      { "F0", required_argument, NULL, 'F' },
	                                      { "fr", required_argument, NULL, 'r' },
	                                      { "f1", required_argument, NULL, 'f' },
	                                      { "prefix", required_argument, NULL, 'p' },
	                                      { "outbloom", no_argument, &opt::outbloom, 1 },
	                                      { "solid", no_argument, &opt::solid, 1 },
	                                      { "eval", no_argument, &opt::eval, 1 },
	                                      { "help", no_argument, NULL, OPT_HELP },
	                                      { "version", no_argument, NULL, OPT_VERSION },
	                                      { NULL, 0, NULL, 0 } };

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

struct entry
{
	string kmer;
	unsigned count;
};

bool
hitSearchInsert(const uint64_t kint, const string& kmer, omp_lock_t* locks, entry* T)
{
	uint64_t i = 0, j;
	do {
		j = (kint + i) % opt::hitSize;
		if (T[j].kmer == kmer) {
#pragma omp atomic
			++T[j].count;
			return true;
		}
		++i;
	} while (i != opt::hitSize && T[j].count != 0);
	if (T[j].count == 0) {
		omp_set_lock(&locks[(uint16_t)j]);
		T[j].kmer = kmer;
		++T[j].count;
		omp_unset_lock(&locks[(uint16_t)j]);
		return false;
	}
	return false;
}

unsigned
hitSearch(const uint64_t kint, const string& kmer, entry* T)
{
	uint64_t i = 0, j;
	do {
		j = (kint + i) % opt::hitSize;
		if (T[j].kmer == kmer) {
			return T[j].count;
		}
		++i;
	} while (i != opt::hitSize && T[j].count != 0);
	return 0;
}

void
fqHit(std::ifstream& in, omp_lock_t* locks, BloomFilter& mydBF, CBFilter& mycBF, entry* hitTable)
{
	bool good, good2 = true;
#pragma omp parallel
	for (string rSeqs, hseq; good2;) {
#pragma omp critical(in)
		{
			good = static_cast<bool>(getline(in, rSeqs));
			good = static_cast<bool>(getline(in, hseq));
			good = static_cast<bool>(getline(in, hseq));
			good2 = static_cast<bool>(getline(in, hseq));
		}
		if (good) {
			ntHashIterator itr(rSeqs, opt::h + 1, opt::k);
			while (itr != itr.end()) {
				if (mydBF.insertAndCheck(*itr)) {
					string canonKmer = rSeqs.substr(itr.pos(), opt::k);
					getCanon(canonKmer);
					if (opt::hitCap > 1) {
						if (mycBF.insert_and_test(*itr)) {
							hitSearchInsert((*itr)[0], canonKmer, locks, hitTable);
						}
					} else {
						hitSearchInsert((*itr)[0], canonKmer, locks, hitTable);
					}
				}
				++itr;
			}
		}
	}
}

void
faHit(std::ifstream& in, omp_lock_t* locks, BloomFilter& mydBF, CBFilter& mycBF, entry* hitTable)
{
	bool good = true;
#pragma omp parallel
	for (string seq, hseq; good;) {
		string rSeqs;
#pragma omp critical(in)
		{
			good = static_cast<bool>(getline(in, seq));
			while (good && seq[0] != '>') {
				rSeqs += seq;
				good = static_cast<bool>(getline(in, seq));
			}
		}
		ntHashIterator itr(rSeqs, opt::h, opt::k);
		while (itr != itr.end()) {
			if (mydBF.insertAndCheck(*itr))
				if (mycBF.insert_and_test(*itr)) {
					string canonKmer = rSeqs.substr(itr.pos(), opt::k);
					getCanon(canonKmer);
					hitSearchInsert((*itr)[0], canonKmer, locks, hitTable);
				}
			++itr;
		}
	}
}

void
bfqHit(std::ifstream& in, BloomFilter& mydBF, CBFilter& mycBF, BloomFilter& myhBF)
{
	bool good, good2 = true;
#pragma omp parallel
	for (string rSeqs, hseq; good2;) {
#pragma omp critical(in)
		{
			good = static_cast<bool>(getline(in, rSeqs));
			good = static_cast<bool>(getline(in, hseq));
			good = static_cast<bool>(getline(in, hseq));
			good2 = static_cast<bool>(getline(in, hseq));
		}
		if (good) {
			ntHashIterator itr(rSeqs, opt::h + 1, opt::k);
			while (itr != itr.end()) {
				if (mydBF.insertAndCheck(*itr)) {
					if (opt::hitCap > 1) {
						if (mycBF.insert_and_test(*itr))
							myhBF.insert(*itr);
					} else {
						myhBF.insert(*itr);
					}
				}
				++itr;
			}
		}
	}
}

void
bfaHit(std::ifstream& in, BloomFilter& mydBF, CBFilter& mycBF, BloomFilter& myhBF)
{
	bool good = true;
#pragma omp parallel
	for (string seq, hseq; good;) {
		string rSeqs;
#pragma omp critical(in)
		{
			good = static_cast<bool>(getline(in, seq));
			while (good && seq[0] != '>') {
				rSeqs += seq;
				good = static_cast<bool>(getline(in, seq));
			}
		}
		ntHashIterator itr(rSeqs, opt::h + 1, opt::k);
		while (itr != itr.end()) {
			if (mydBF.insertAndCheck(*itr)) {
				if (opt::hitCap > 1) {
					if (mycBF.insert_and_test(*itr))
						myhBF.insert(*itr);
				} else {
					myhBF.insert(*itr);
				}
			}
			++itr;
		}
	}
}

int
main(int argc, char** argv)
{
	double sTime = omp_get_wtime();

	bool die = false;
	bool nontCard = false;
	for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
		std::istringstream arg(optarg != NULL ? optarg : "");
		switch (c) {
		case '?':
			die = true;
			break;
		case 't':
			arg >> opt::t;
			break;
		case 'h':
			arg >> opt::h;
			break;
		case 'k':
			arg >> opt::k;
			break;
		case 'b':
			arg >> opt::m;
			break;
		case 'c':
			arg >> opt::hitCap;
			break;
		case 'F':
			arg >> opt::F0;
			nontCard = true;
			break;
		case 'r':
			arg >> opt::fr;
			nontCard = true;
			break;
		case 'f':
			arg >> opt::f1;
			nontCard = true;
			break;
		case 'p':
			arg >> opt::prefix;
			break;
		case OPT_HELP:
			std::cerr << USAGE_MESSAGE;
			exit(EXIT_SUCCESS);
		case OPT_VERSION:
			std::cerr << VERSION_MESSAGE;
			exit(EXIT_SUCCESS);
		}
		if (optarg != NULL && !arg.eof()) {
			std::cerr << PROGRAM ": invalid option: `-" << (char)c << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}
	if (argc - optind < 1) {
		std::cerr << PROGRAM ": missing arguments\n";
		die = true;
	}

	if (opt::k == 0) {
		std::cerr << PROGRAM ": missing argument -k ... \n";
		die = true;
	}

	if (die) {
		std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	vector<string> inFiles;
	for (int i = optind; i < argc; ++i) {
		string file(argv[i]);
		if (file[0] == '@') {
			string inName;
			ifstream inList(file.substr(1, file.length()).c_str());
			while (getline(inList, inName))
				inFiles.push_back(inName);
		} else
			inFiles.push_back(file);
	}
	if (!nontCard) {
		size_t histArray[10002];
		getHist(inFiles, opt::k, opt::t, histArray);

		int histIndex = 2, errCov = 1;
		while (histIndex <= 10000 && histArray[histIndex] > histArray[histIndex + 1])
			histIndex++;

		errCov = histIndex > 300 ? 1 : histIndex - 1;

		unsigned maxCov = errCov;
		for (unsigned i = errCov; i < 10002; i++) {
			if (histArray[i] >= histArray[maxCov])
				maxCov = i;
		}
		maxCov--;

		unsigned minCov = errCov;
		for (unsigned i = errCov; i < maxCov; i++) {
			if (histArray[i] <= histArray[minCov])
				minCov = i;
		}
		minCov--;

		if (opt::solid) {
			if (opt::hitCap == 0)
				opt::hitCap = minCov;
			cerr << "Errors k-mer coverage: " << opt::hitCap << endl;
		}

		if (!opt::solid) {
			if (opt::hitCap == 0)
				opt::hitCap = 1.75 * maxCov;
			cerr << "Errors k-mer coverage: " << minCov << endl;
			cerr << "Median k-mer coverage: " << maxCov << endl;
			cerr << "Repeat k-mer coverage: " << opt::hitCap << endl;
		}

		opt::dbfSize = opt::bits * histArray[1];
		opt::cbfSize = opt::bytes * (histArray[1] - histArray[2]);
		size_t hitCount = histArray[1];
		for (unsigned i = 2; i <= opt::hitCap + 1; i++)
			hitCount -= histArray[i];
		opt::hitSize = opt::outbloom ? hitCount * opt::m : hitCount * 3;

		cerr << "Approximate# of distinct k-mers: " << histArray[1] << "\n";
		cerr << "Approximate# of solid k-mers: " << hitCount << "\n";
	} else {
		opt::dbfSize = opt::bits * opt::F0;
		opt::cbfSize = opt::bytes * (opt::F0 - opt::f1);
		opt::hitSize = opt::m * opt::fr;
		cerr << "Approximate# of distinct k-mers: " << opt::F0 << "\n";
		cerr << "Approximate# of solid k-mers: " << opt::fr << "\n";
	}

#ifdef _OPENMP
	omp_set_num_threads(opt::t);
#endif

	BloomFilter mydBF(((opt::dbfSize + 7) / 8) * 8, 3, opt::k);
	std::cerr << "mydBF size: " << mydBF.getFilterSize() << std::endl;
	CBFilter mycBF(opt::cbfSize, opt::h, opt::k, opt::hitCap - 1);

	if (opt::outbloom) {
		BloomFilter myhBF(((opt::hitSize + 7) / 8) * 8, opt::h + 1, opt::k);
		std::cerr << "myhBF size: " << myhBF.getFilterSize() << std::endl;
		for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
			std::ifstream in(inFiles[file_i].c_str());
			string firstLine;
			bool good = static_cast<bool>(getline(in, firstLine));
			if (!good) {
				std::cerr << "Error in reading file: " << inFiles[file_i] << "\n";
				exit(EXIT_FAILURE);
			}
			if (firstLine[0] == '>')
				bfaHit(in, mydBF, mycBF, myhBF);
			else if (firstLine[0] == '@')
				bfqHit(in, mydBF, mycBF, myhBF);
			else {
				std::cerr << "Error in reading file: " << inFiles[file_i] << "\n";
				exit(EXIT_FAILURE);
			}
			in.close();
		}

		std::stringstream hstm;
		if (opt::prefix.empty()) {
			if (opt::solid)
				hstm << "solids_k" << opt::k << ".bf";
			else
				hstm << "repeat_k" << opt::k << ".bf";
		} else
			hstm << opt::prefix << "_k" << opt::k << ".bf";
		myhBF.storeFilter(hstm.str().c_str());
	} else {
		entry* hitTable = new entry[opt::hitSize];
		for (size_t i = 0; i < opt::hitSize; i++)
			hitTable[i].count = 0;
		const unsigned lockSize = 65536;
		omp_lock_t* locks = new omp_lock_t[lockSize];

		for (unsigned i = 0; i < lockSize; i++)
			omp_init_lock(&locks[i]);

		for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
			std::ifstream in(inFiles[file_i].c_str());
			string firstLine;
			bool good = static_cast<bool>(getline(in, firstLine));
			if (!good) {
				std::cerr << "Error in reading file: " << inFiles[file_i] << "\n";
				exit(EXIT_FAILURE);
			}
			if (firstLine[0] == '>')
				faHit(in, locks, mydBF, mycBF, hitTable);
			else if (firstLine[0] == '@')
				fqHit(in, locks, mydBF, mycBF, hitTable);
			else {
				std::cerr << "Error in reading file: " << inFiles[file_i] << "\n";
				exit(EXIT_FAILURE);
			}
			in.close();
		}

		for (unsigned i = 0; i < lockSize; i++)
			omp_destroy_lock(&locks[i]);
		delete[] locks;

		std::stringstream hstm;
		if (opt::prefix.empty()) {
			if (opt::solid)
				hstm << "solids_k" << opt::k << ".rep";
			else
				hstm << "repeat_k" << opt::k << ".rep";
		} else
			hstm << opt::prefix << "_k" << opt::k << ".rep";
		ofstream outFile(hstm.str().c_str());
		for (size_t i = 0; i < opt::hitSize; i++)
			if (hitTable[i].count != 0)
				outFile << hitTable[i].kmer << "\t" << hitTable[i].count + opt::hitCap << "\n";
		outFile.close();

		// if(opt::eval) {
		std::ofstream dOut("diff_2");
		std::ifstream dskRes("out_ge22.txt");
		std::string dskLine;
		while (getline(dskRes, dskLine)) {
			std::string dsk_kmer;
			unsigned dsk_count = 0;
			std::istringstream dskstm(dskLine);
			dskstm >> dsk_kmer >> dsk_count;
			// cerr << dsk_kmer << "\t" << dsk_count << "\n";
			// if (dsk_count >= opt::hitCap) {
			uint64_t dsk_hash = NTC64(dsk_kmer.c_str(), opt::k);
			unsigned nthits_count = hitSearch(dsk_hash, dsk_kmer, hitTable);
			if (nthits_count != 0) {
				std::cout << dsk_kmer << "\t" << dsk_count << "\t" << nthits_count + opt::hitCap
				          << "\n";
			} else {
				std::string dsk_can(dsk_kmer);
				getCanon(dsk_can);
				if (dsk_can != dsk_kmer) {
					uint64_t dsk_hash = NTC64(dsk_can.c_str(), opt::k);
					unsigned nthits_count = hitSearch(dsk_hash, dsk_can, hitTable);
					if (nthits_count != 0) {
						std::cout << dsk_can << "\t" << dsk_count << "\t"
						          << nthits_count + opt::hitCap << "\n";
					}
				} else
					dOut << dsk_kmer << "\t" << dsk_count << "\t"
					     << "0"
					     << "\n";
			}
			//}
		}

		dskRes.close();
		dOut.close();
		//}

		delete[] hitTable;
	}

	cerr << "Total time for computing repeat content in (sec): " << setprecision(4) << fixed
	     << omp_get_wtime() - sTime << "\n";
	return 0;
}
