#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "ntHashIterator.hpp"
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
    unsigned nThrd=8;
    unsigned k = 32;
    size_t cbfSize = 110000000;
    size_t hitSize = 1000000;
    unsigned hitCap = 20;
    size_t counter = 0;
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
            T[j].count = T[j].count+1;
            return true;
        }
        ++i;
    } while (i!= opt::hitSize && T[j].count!=0);
    if (T[i].count == 0) {
        #pragma omp critical (hashUpdate)
        {
            T[j].kmer = kmer;
            T[j].count = 1;
        }
        return false;
    }
    return false;
}

int main(int argc, char** argv) {

    double sTime = omp_get_wtime();

    opt::hitSize = 36000000;
    opt::cbfSize = 3100000000;
    opt::hitCap = 40;
    entry *hitTable = new entry [opt::hitSize];
    for (unsigned i=0; i<opt::hitSize; i++) {
        //hitTable[i].kmer = "";
        hitTable[i].count = 0;
    }

    CBFilter myCBF(opt::cbfSize, 3, opt::k, opt::hitCap);
    ifstream in(argv[1]);
    
#ifdef _OPENMP
    omp_set_num_threads(opt::nThrd);
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
            ntHashIterator itr(rSeqs, 3, opt::k);
            while (itr != itr.end()) {
                if(myCBF.insert_and_test(*itr))
                    hitSearchInsert((*itr)[0], rSeqs.substr(itr.pos(),opt::k), hitTable);
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
