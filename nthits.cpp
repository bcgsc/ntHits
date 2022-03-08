#include <cmath>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "CBFilter.hpp"
#include "btllib/bloom_filter.hpp"
#include "btllib/nthash.hpp"
#include "btllib/seq_reader.hpp"
#include "ntcard.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "nthits"

static const char VERSION_MESSAGE[] =
  PROGRAM " 0.1.0 \n"
          "Written by Hamid Mohamadi.\n"
          "Copyright 2019 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
  "Usage: " PROGRAM " [OPTION]... FILES...\n"
  "Reports the most frequent k-mers in FILES(>=1).\n"
  "Acceptable file formats: fastq, fasta, gz, bz, zip.\n"
  "\n"
  " Options:\n"
  "\n"
  "  -t, --threads=N	use N parallel threads [16]\n"
  "  -k, --kmer=N	the length of kmer [64]\n"
  "  -c, --cutoff=N	the maximum coverage of kmer in output\n"
  "  -p, --pref=STRING	the prefix for output file name [repeat]\n"
  "  --outbloom	output the most frequent k-mers in a Bloom filter\n"
  "  --solid	output the solid k-mers (non-erroneous k-mers)\n"
	"  --long   enable long mode for optimized reading of long reads\n"
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
int longmode;
} // namespace opt

static const char shortopts[] = "t:h:k:b:p:r:c:F:f:";

enum
{
  OPT_HELP = 1,
  OPT_VERSION
};

static const struct option longopts[] = {
  { "threads", required_argument, NULL, 't' },
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
  { "long", no_argument, &opt::longmode, 1 },
  { "help", no_argument, NULL, OPT_HELP },
  { "version", no_argument, NULL, OPT_VERSION },
  { NULL, 0, NULL, 0 }
};

static const unsigned char b2r[256] = {
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', // 0
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
  'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

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
    for (int lIndex = p, rIndex = bMer.length() - 1 - p; lIndex <= rIndex;
         ++lIndex, --rIndex) {
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
hitSearchInsert(const uint64_t kint,
                const string& kmer,
                omp_lock_t* locks,
                entry* T)
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
fHit(const std::string& filepath,
     omp_lock_t* locks,
     btllib::KmerBloomFilter& mydBF,
     CBFilter& mycBF,
     entry* hitTable)
{
  const auto flags = opt::longmode ? btllib::SeqReader::Flag::LONG_MODE
                                   : btllib::SeqReader::Flag::SHORT_MODE;
  btllib::SeqReader reader(filepath, flags);
#pragma omp parallel shared(reader)
  for (const auto record : reader) {
    btllib::NtHash nthash(record.seq, opt::h + 1, opt::k);
    while (nthash.roll()) {
      if (mydBF.contains_insert(nthash.hashes())) {
        std::string canonKmer = record.seq.substr(nthash.get_pos(), opt::k);
        getCanon(canonKmer);
        if (opt::hitCap > 1) {
          if (mycBF.insert_and_test(nthash.hashes())) {
            hitSearchInsert(nthash.hashes()[0], canonKmer, locks, hitTable);
          }
        } else {
          hitSearchInsert(nthash.hashes()[0], canonKmer, locks, hitTable);
        }
      }
    }
  }
}

void
bHit(const std::string& filepath,
     btllib::KmerBloomFilter& mydBF,
     CBFilter& mycBF,
     btllib::KmerBloomFilter& myhBF)
{
  const auto flags = opt::longmode ? btllib::SeqReader::Flag::LONG_MODE
                                   : btllib::SeqReader::Flag::SHORT_MODE;
  btllib::SeqReader reader(filepath, flags);
#pragma omp parallel shared(reader)
  for (const auto record : reader) {
    btllib::NtHash nthash(record.seq, opt::h + 1, opt::k);
    while (nthash.roll()) {
      if (mydBF.contains_insert(nthash.hashes())) {
        if (opt::hitCap > 1) {
          if (mycBF.insert_and_test(nthash.hashes()))
            myhBF.insert(nthash.hashes());
        } else {
          myhBF.insert(nthash.hashes());
        }
      }
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
    while (histIndex <= 10000 &&
           histArray[histIndex] > histArray[histIndex + 1])
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
  CBFilter mycBF(opt::cbfSize, opt::h, opt::k, opt::hitCap - 1);

  if (opt::outbloom) {
    btllib::KmerBloomFilter myhBF(opt::hitSize / 8, opt::h + 1, opt::k);
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
      bHit(inFiles[file_i], mydBF, mycBF, myhBF);
    }

    std::stringstream hstm;
    if (opt::prefix.empty()) {
      if (opt::solid)
        hstm << "solids_k" << opt::k << ".bf";
      else
        hstm << "repeat_k" << opt::k << ".bf";
    } else
      hstm << opt::prefix << "_k" << opt::k << ".bf";
    myhBF.save(hstm.str());
  } else {
    entry* hitTable = new entry[opt::hitSize];
    for (size_t i = 0; i < opt::hitSize; i++)
      hitTable[i].count = 0;
    const unsigned lockSize = 65536;
    omp_lock_t* locks = new omp_lock_t[lockSize];

    for (unsigned i = 0; i < lockSize; i++)
      omp_init_lock(&locks[i]);

    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
      fHit(inFiles[file_i], locks, mydBF, mycBF, hitTable);
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
        outFile << hitTable[i].kmer << "\t" << hitTable[i].count + opt::hitCap
                << "\n";
    outFile.close();

    delete[] hitTable;
  }

  cerr << "Total time for computing repeat content in (sec): "
       << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
  return 0;
}
