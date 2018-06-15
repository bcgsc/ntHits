#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "ntHashIterator.hpp"
#include "CBFilter.hpp"


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
    unsigned k = 32;
    uint64_t relZero = 0;
    size_t hitSize = 0;
    size_t maxSize = 1000000;
    size_t counter=0;
}



struct entry {
    string kmer;
    size_t count;
};

bool hitSearchInsert(const uint64_t kint, const string &kmer, entry *T){
    ++opt::counter;
    uint64_t i=0, j;
    do {
        j = (kint + i) % opt::maxSize;
        if (T[j].kmer == kmer) {
            T[j].count = T[j].count+1;
            //cerr << opt::counter << "\t" << i << ":" << j << "\ttrue\t" << T[j].kmer << "\t" << T[j].count << endl;
            return true;
        }
        ++i;
    }while (i!= opt::maxSize && T[j].count!=0);
    if (T[i].count == 0) {
        T[j].kmer = kmer;
        T[j].count = 1;
        return false;
    }
    return false;
}

int main(int argc, char** argv) {

    clock_t sTime = clock();

    //size_t opt::maxSize = ntCard(file);
    size_t F0 = 110000000;
    entry *hitTable = new entry [opt::maxSize];
    for (unsigned i=0; i<opt::maxSize; i++) {
        //hitTable[i].kmer = "";
        hitTable[i].count = 0;
    }

    CBFilter myCBF(F0, 3, opt::k); 
    ifstream mFile(argv[1]);
    
    string rHead,rSeqs, rDirc, rQual;
    while(getline(mFile,rHead)&&getline(mFile,rSeqs)&&getline(mFile,rDirc)&&getline(mFile, rQual)) {
        ntHashIterator itr(rSeqs, 3, opt::k);
        while (itr != itr.end()) {
            //string myKmer = rSeqs.substr(i,opt::k);
            if(myCBF.insert_and_test(*itr))
                hitSearchInsert((*itr)[0], rSeqs.substr(itr.pos(),opt::k), hitTable);
            ++itr;
        }
    }
    mFile.close();
    
    ofstream outFile("hmers");
    for (unsigned i=0; i<opt::maxSize; i++)
        if(hitTable[i].count != 0)
            outFile << hitTable[i].kmer << "\t" << hitTable[i].count << "\n";
    
    //cerr << "relZero=" << opt::relZero << endl;

    
    outFile.close();
    delete [] hitTable;
    
    std::cerr << (double)(clock() - sTime)/CLOCKS_PER_SEC << "\n";
    return 0;
}
