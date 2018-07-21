#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <getopt.h>


#include "ntHashIterator.hpp"
#include "ntcard.hpp"
#include "CBFilter.hpp"

#include "Uncompress.h"

#ifdef _OPENMP
# include <omp.h>
#endif


#define PROGRAM "nthits"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 0.0.1 \n"
    "Written by Hamid Mohamadi.\n"
    "Copyright 2018 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... FILES...\n"
    "Reports the most frequent k-mers in FILES(>=1).\n"
    "Accepatble file formats: fastq, fasta, gz, bz, zip.\n"
    "\n"
    " Options:\n"
    "\n"
    "  -j, --threads=N	use N parallel threads [16]\n"
    "  -k, --kmer=N	the length of kmer [64]\n"
    "  -c, --cutoff=N	the maximum coverage of kmer in output [40]\n"
    "  -b, --bit=N	the number of counter per distict k-mer [64]\n"
    "  -p, --pref=STRING	the prefix for output file name [repeat]\n"

    "      --help	display this help and exit\n"
    "      --version	output version information and exit\n"
    "\n"
    "Report bugs to hmohamadi@bcgsc.ca.\n";


using namespace std;

namespace opt {
unsigned j = 16;
unsigned k = 64;
unsigned h = 2;
unsigned hitCap=40;
unsigned bits = 3;
size_t cbfSize;
size_t hitSize;
string prefix;
}

static const char shortopts[] = "j:h:k:b::p:c:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 'j' },
    { "kmer",	required_argument, NULL, 'k' },
    { "hash",	required_argument, NULL, 'h' },
    { "bit",	required_argument, NULL, 'b' },
    { "cutoff",	required_argument, NULL, 'c' },
    { "prefix",	required_argument, NULL, 'p' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};


static const unsigned char b2r[256] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //0
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //1
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //2
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //3
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', //4   'A' 'C' 'G'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', //5   'T'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', //6   'a' 'c' 'g'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', //7   't'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //8
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //9
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //10
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //11
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //12
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //13
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //14
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //15
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

void getCanon(std::string &bMer) {
    int p=0, hLen=(bMer.length()-1)/2;
    while (bMer[p] == b2r[(unsigned char)bMer[bMer.length()-1-p]]) {
        ++p;
        if(p>=hLen) break;
    }
    if (bMer[p] > b2r[(unsigned char)bMer[bMer.length()-1-p]]) {
        for (int lIndex = p, rIndex = bMer.length()-1-p; lIndex<=rIndex; ++lIndex,--rIndex) {
            char tmp = b2r[(unsigned char)bMer[rIndex]];
            bMer[rIndex] = b2r[(unsigned char)bMer[lIndex]];
            bMer[lIndex] = tmp;
        }
    }
}

struct entry {
    string kmer;
    size_t count;
};

bool hitSearchInsert(const uint64_t kint, const string &kmer, omp_lock_t *locks, entry *T) {
    uint64_t i=0, j;
    do {
        j = (kint + i) % opt::hitSize;
        if (T[j].kmer == kmer) {
            #pragma omp atomic
            ++T[j].count;
            return true;
        }
        ++i;
    } while (i!= opt::hitSize && T[j].count!=0);
    if (T[j].count == 0) {
        omp_set_lock(&locks[(uint16_t)j]);
        T[j].kmer = kmer;
        ++T[j].count;
        omp_unset_lock(&locks[(uint16_t)j]);
        return false;
    }
    return false;
}

void fqHit(std::ifstream &in, omp_lock_t *locks, CBFilter &myCBF, entry *hitTable) {
    bool good, good2 =true;
    #pragma omp parallel
    for(string rSeqs, hseq; good2;) {
        #pragma omp critical(in)
        {
            good = static_cast<bool>(getline(in, rSeqs));
            good = static_cast<bool>(getline(in, hseq));
            good = static_cast<bool>(getline(in, hseq));
            good2 = static_cast<bool>(getline(in, hseq));
        }
        if(good) {
            ntHashIterator itr(rSeqs, opt::h, opt::k);
            while (itr != itr.end()) {
                if(myCBF.insert_and_test(*itr)) {
                    string canonKmer = rSeqs.substr(itr.pos(),opt::k);
                    getCanon(canonKmer);
                    hitSearchInsert((*itr)[0], canonKmer, locks, hitTable);
                }
                ++itr;
            }
        }
    }
}

void faHit(std::ifstream &in, omp_lock_t *locks, CBFilter &myCBF, entry *hitTable) {
    bool good = true;
    #pragma omp parallel
    for(string seq, hseq; good;) {
        string rSeqs;
        #pragma omp critical(in)
        {
            good = static_cast<bool>(getline(in, seq));
            while(good&&seq[0]!='>') {
                rSeqs+=seq;
                good = static_cast<bool>(getline(in, seq));
            }
        }
        ntHashIterator itr(rSeqs, opt::h, opt::k);
        while (itr != itr.end()) {
            if(myCBF.insert_and_test(*itr)) {
                string canonKmer = rSeqs.substr(itr.pos(),opt::k);
                getCanon(canonKmer);
                hitSearchInsert((*itr)[0], canonKmer, locks, hitTable);
            }
            ++itr;
        }
    }
}

int main(int argc, char** argv) {

    double sTime = omp_get_wtime();

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 'j':
            arg >> opt::j;
            break;
        case 'h':
            arg >> opt::h;
            break;
        case 'k':
            arg >> opt::k;
            break;
        case 'b':
            arg >> opt::bits;
            break;
        case 'c':
            arg >> opt::hitCap;
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
            std::cerr << PROGRAM ": invalid option: `-"
                      << (char)c << optarg << "'\n";
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
        if(file[0]=='@') {
            string inName;
            ifstream inList(file.substr(1,file.length()).c_str());
            while(getline(inList,inName))
                inFiles.push_back(inName);
        }
        else
            inFiles.push_back(file);
    }

    size_t histArray[10002];
    getHist(inFiles, opt::k, opt::j, histArray);

    opt::cbfSize = opt::bits*histArray[1];
    opt::hitSize = histArray[1];
    for(unsigned i=2; i<=opt::hitCap; i++)
        opt::hitSize -= histArray[i];
    opt::hitSize *= 5;


    entry *hitTable = new entry [opt::hitSize];
    for (unsigned i=0; i<opt::hitSize; i++) {
        hitTable[i].kmer = "";
        hitTable[i].count = 0;
    }

    CBFilter myCBF(opt::cbfSize, opt::h, opt::k, opt::hitCap);

#ifdef _OPENMP
    omp_set_num_threads(opt::j);
#endif

    unsigned lockSize = 65536;
    omp_lock_t *locks = new omp_lock_t [lockSize];
    for(unsigned i = 0; i < lockSize; i++) omp_init_lock(&locks[i]);

    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
        std::ifstream in(inFiles[file_i].c_str());
        string firstLine;
        bool good = static_cast<bool>(getline(in, firstLine));
        if (!good) {
            std::cerr << "Error in reading file: " << inFiles[file_i] << "\n";
            exit(EXIT_FAILURE);
        }
        if(firstLine[0]=='>')
            faHit(in, locks, myCBF, hitTable);
        else if (firstLine[0]=='@')
            fqHit(in, locks, myCBF, hitTable);
        else {
            std::cerr << "Error in reading file: " << inFiles[file_i] << "\n";
            exit(EXIT_FAILURE);
        }
        in.close();
    }

    for(unsigned i = 0; i < lockSize; i++) omp_destroy_lock(&locks[i]);
    delete [] locks;


    std::stringstream hstm;
    if(opt::prefix.empty())
        hstm << "repeat_k" << opt::k << ".rep";
    else
        hstm << opt::prefix << "_k" << opt::k << ".rep";
    ofstream outFile(hstm.str().c_str());
    for (unsigned i=0; i<opt::hitSize; i++)
        if(hitTable[i].count != 0)
            outFile << hitTable[i].kmer << "\t" << hitTable[i].count + opt::hitCap  << "\n";
    outFile.close();

    delete [] hitTable;

    cerr << "Total time for computing reapeat content in (sec): " <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    return 0;
}
