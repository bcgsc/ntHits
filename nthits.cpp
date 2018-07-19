#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "ntHashIterator.hpp"
#include "ntcard.hpp"
#include "CBFilter.hpp"

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
    "Estimates the number of k-mers in FILES(>=1).\n"
    "Accepatble file formats: fastq, fasta, sam, bam, gz, bz, zip.\n"
    "\n"
    " Options:\n"
    "\n"
    "      --help	display this help and exit\n"
    "      --version	output version information and exit\n"
    "\n"
    "Report bugs to hmohamadi@bcgsc.ca.\n";


using namespace std;

namespace opt {
    unsigned j = 24;
    unsigned k = 64;
    unsigned h = 3;
    unsigned bits = 5;
    size_t cbfSize;
    size_t hitSize;
    unsigned hitCap;
}


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

bool hitSearchInsert(const uint64_t kint, const string &kmer, entry *T){
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
        #pragma omp critical (T)
        {
            T[j].kmer = kmer;
            ++T[j].count;
        }
        return false;
    }
    return false;
}

int main(int argc, char** argv) {

    double sTime = omp_get_wtime();
    
    vector<string> inFiles;
    inFiles.push_back(argv[1]);
    
    opt::h = atoi(argv[2]);
    opt::bits = atoi(argv[3]);
    size_t histArray[10002];
    getHist(inFiles, opt::k, opt::j, histArray);
    
    opt::cbfSize = opt::bits*histArray[1];//3*histArray[1]; // *5=
    opt::hitCap = 40;
    //opt::hitSize = 25522210;
    opt::hitSize = histArray[1];
    for(unsigned i=2; i<=opt::hitCap; i++)
        opt::hitSize -= histArray[i];
    opt::hitSize *= 5;
    
	cerr << opt::cbfSize << "\t" << opt::hitSize << "\n\n\n";	

    entry *hitTable = new entry [opt::hitSize];
    for (unsigned i=0; i<opt::hitSize; i++) {
        hitTable[i].count = 0;
    }

    CBFilter myCBF(opt::cbfSize, opt::h, opt::k, opt::hitCap);
    ifstream in(argv[1]);
    
#ifdef _OPENMP
    omp_set_num_threads(opt::j);
#endif

    bool good = true;
    #pragma omp parallel
    for(string rSeqs, hseq; good;) {
        #pragma omp critical(in)
        {
            good = static_cast<bool>(getline(in, hseq));
            good = static_cast<bool>(getline(in, rSeqs));
            good = static_cast<bool>(getline(in, hseq));
            good = static_cast<bool>(getline(in, hseq));
        }
        if(good) {
            ntHashIterator itr(rSeqs, opt::h, opt::k);
            while (itr != itr.end()) {
                if(myCBF.insert_and_test(*itr)) {
                    string canonKmer = rSeqs.substr(itr.pos(),opt::k);
                    getCanon(canonKmer);
                    hitSearchInsert((*itr)[0], canonKmer, hitTable);
                }
                ++itr;
            }
        }
    }
    in.close();
    
    ofstream outFile("hmers");
    for (unsigned i=0; i<opt::hitSize; i++)
        if(hitTable[i].count != 0)
            outFile << hitTable[i].kmer << "\t" << hitTable[i].count + opt::hitCap  << "\n";
    
    outFile.close();
    delete [] hitTable;
    
    cerr << "time(sec): " <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    return 0;
}
