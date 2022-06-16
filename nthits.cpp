#include <argparse/argparse.hpp>
#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <btllib/nthash.hpp>
#include <btllib/seq_reader.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "ntcard.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM_NAME "ntHits"
#define PROGRAM_VERSION "0.0.1"
#define PROGRAM_DESCRIPTION "Reports the most frequent k-mers in input files."
#define PROGRAM_COPYRIGHT "Copyright 2019 Canada's Michael Smith Genome Science Centre"

using namespace std;
using cbf_counter_t = uint8_t;

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
bool outbloom = false;
bool solid = false;
bool eval = false;
unsigned seq_reader_mode;
} // namespace opt

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
f_hit(
    const std::string& file_path,
    omp_lock_t* locks,
    btllib::KmerBloomFilter& mydBF,
    btllib::CountingBloomFilter<cbf_counter_t>& mycBF,
    entry* hitTable)
{
	btllib::SeqReader reader(file_path, opt::seq_reader_mode);
#pragma omp parallel shared(reader)
	for (const auto& record : reader) {
		btllib::NtHash nth(record.seq, opt::h + 1, opt::k);
		while (nth.roll()) {
			if (!mydBF.contains_insert(nth.hashes())) {
				string canonKmer = record.seq.substr(nth.get_pos(), opt::k);
				getCanon(canonKmer);
				if (opt::hitCap > 1) {
					if (mycBF.insert_thresh_contains(nth.hashes(), opt::hitCap - 1)) {
						hitSearchInsert(nth.hashes()[0], canonKmer, locks, hitTable);
					}
				} else {
					hitSearchInsert(nth.hashes()[0], canonKmer, locks, hitTable);
				}
			}
		}
	}
}

void
b_hit(
    const std::string& file_path,
    btllib::KmerBloomFilter& mydBF,
    btllib::CountingBloomFilter<cbf_counter_t>& mycBF,
    btllib::KmerBloomFilter& myhBF)
{
	btllib::SeqReader reader(file_path, opt::seq_reader_mode);
#pragma omp parallel shared(reader)
	for (const auto record : reader) {
		btllib::NtHash nth(record.seq, opt::h + 1, opt::k);
		while (nth.roll()) {
			if (mydBF.contains_insert(nth.hashes())) {
				if (opt::hitCap > 1) {
					if (mycBF.insert_thresh_contains(nth.hashes(), opt::hitCap - 1))
						myhBF.insert(nth.hashes());
				} else {
					myhBF.insert(nth.hashes());
				}
			}
		}
	}
}

int
main(int argc, char** argv)
{
	double sTime = omp_get_wtime();

	bool nontCard = false;
	std::vector<std::string> inFiles;

	auto default_args = argparse::default_arguments::none;
	argparse::ArgumentParser parser(PROGRAM_NAME, PROGRAM_VERSION, default_args);
	parser.add_description(PROGRAM_DESCRIPTION);
	parser.add_epilog(PROGRAM_COPYRIGHT);

	parser.add_argument("-t", "--threads")
	    .help("Number of parallel threads")
	    .default_value(16U)
	    .scan<'u', unsigned>();

	parser.add_argument("-k", "--kmer")
	    .help("k-mer length")
	    .default_value(64U)
	    .scan<'u', unsigned>();

	parser.add_argument("-h", "--hashes")
	    .help("Number of hashes to generate per k-mer/spaced seed")
	    .default_value(4U)
	    .scan<'u', unsigned>();

	parser.add_argument("-c", "--cutoff")
	    .help("k-mer cutoff threshold")
	    .required()
	    .scan<'u', unsigned>();

	parser.add_argument("-p", "--prefix")
	    .help("Output files' prefix")
	    .default_value(std::string("repeat"));

	parser.add_argument("--outbloom")
	    .help("Output the most frequent k-mers in a Bloom filter")
	    .default_value(false)
	    .implicit_value(true);

	parser.add_argument("--solid")
	    .help("Output the solid k-mers (non-erroneous k-mers)")
	    .default_value(false)
	    .implicit_value(true);

	parser.add_argument("--long-mode")
	    .help("Optimize data reader for long sequences (>5kbp)")
	    .default_value(false)
	    .implicit_value(true);

	parser.add_argument("-b", "--bit").default_value(16U).scan<'u', unsigned>();
	parser.add_argument("-F").scan<'u', unsigned>();
	parser.add_argument("-f").scan<'u', unsigned>();
	parser.add_argument("-r").scan<'u', unsigned>();

	parser.add_argument("files").help("Input files").required().remaining();

	try {
		parser.parse_args(argc, argv);
	} catch (const std::runtime_error& err) {
		std::cerr << err.what() << std::endl;
		std::cerr << parser;
		std::exit(1);
	}

	opt::t = parser.get<unsigned>("-t");
	opt::k = parser.get<unsigned>("-k");
	opt::h = parser.get<unsigned>("-h");
	opt::m = parser.get<unsigned>("-b");
	opt::hitCap = parser.get<unsigned>("-c");
	opt::prefix = parser.get("-p");
	opt::outbloom = parser.get<bool>("--outbloom");
	opt::solid = parser.get<bool>("--solid");

	if (parser.get<bool>("--long-mode")) {
		opt::seq_reader_mode = btllib::SeqReader::Flag::LONG_MODE;
	} else {
		opt::seq_reader_mode = btllib::SeqReader::Flag::SHORT_MODE;
	}

	if (parser.is_used("-F")) {
		opt::F0 = parser.get<unsigned>("-F");
		nontCard = true;
	}

	if (parser.is_used("-r")) {
		opt::fr = parser.get<unsigned>("-r");
		nontCard = true;
	}

	if (parser.is_used("-f")) {
		opt::f1 = parser.get<unsigned>("-f");
		nontCard = true;
	}

	try {
		inFiles = parser.get<std::vector<std::string> >("files");
	} catch (std::logic_error& e) {
		std::cerr << "No files provided" << std::endl;
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

	btllib::KmerBloomFilter mydBF(opt::dbfSize / 8, 3, opt::k);
	btllib::CountingBloomFilter<cbf_counter_t> mycBF(opt::cbfSize, opt::h);

	if (opt::outbloom) {
		btllib::KmerBloomFilter myhBF(opt::hitSize / 8, opt::h + 1, opt::k);
		for (const auto& file_path : inFiles) {
			b_hit(file_path, mydBF, mycBF, myhBF);
		}
		std::string hbf_out_path;
		if (opt::prefix.empty() && opt::solid) {
			hbf_out_path = "solids_k" + std::to_string(opt::k) + ".bf";
		} else if (opt::prefix.empty() && !opt::solid) {
			hbf_out_path = "repeat_k" + std::to_string(opt::k) + ".bf";
		} else {
			hbf_out_path = opt::prefix + "_k" + std::to_string(opt::k) + ".bf";
		}
		myhBF.save(hbf_out_path);
	} else {
		entry* hitTable = new entry[opt::hitSize];
		for (size_t i = 0; i < opt::hitSize; i++)
			hitTable[i].count = 0;
		const unsigned lockSize = 65536;
		omp_lock_t* locks = new omp_lock_t[lockSize];

		for (unsigned i = 0; i < lockSize; i++)
			omp_init_lock(&locks[i]);

		for (const auto& file_path : inFiles) {
			f_hit(file_path, locks, mydBF, mycBF, hitTable);
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

		delete[] hitTable;
	}

	cerr << "Total time for computing repeat content in (sec): " << setprecision(4) << fixed
	     << omp_get_wtime() - sTime << "\n";
	return 0;
}
