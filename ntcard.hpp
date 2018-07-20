/*
 *
 * ntcard.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */


#ifndef NTCARD_H_
#define NTCARD_H_

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <getopt.h>
#include <cassert>
#include <cmath>
#include "ntHashIterator.hpp"


#ifdef _OPENMP
# include <omp.h>
#endif


using namespace std;

namespace opt {
unsigned nThrd=1;
unsigned kmLen=64;
size_t rBuck;
unsigned rBits=27;
unsigned sBits=11;
unsigned sMask;
unsigned covMax=10000;
size_t nSamp=2;
bool samH=true;
}

static const char shortopts[] = "t:s:r:k:c:l:p:f:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 't' },
    { "kmer",	required_argument, NULL, 'k' },
    { "cov",	required_argument, NULL, 'c' },
    { "rbit",	required_argument, NULL, 'r' },
    { "sbit",	required_argument, NULL, 's' },
    { "pref",	required_argument, NULL, 'p' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

size_t getInf(const char* inFile) {
    std::ifstream in(inFile, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

unsigned getftype(std::ifstream &in, std::string &samSeq) {
    std::string hseq;
    getline(in,hseq);
    if(hseq[0]=='>') {
        return 1;
    }
    if(hseq[0]=='@') {
        if( (hseq[1]=='H'&& hseq[2]=='D') ||
                (hseq[1]=='S'&& hseq[2]=='Q') ||
                (hseq[1]=='R'&& hseq[2]=='G') ||
                (hseq[1]=='P'&& hseq[2]=='G') ||
                (hseq[1]=='C'&& hseq[2]=='O') ) {
            return 2;
        }
        else
            return 0;
    }
    std::istringstream alnSec(hseq);
    std::string s1,s2,s3,s4,s5,s6,s7,s8,s9,s10, s11;
    alnSec>>s1>>s2>>s3>>s4>>s5>>s6>>s7>>s8>>s9>>s10>>s11;
    if((s2.find_first_not_of("0123456789") == std::string::npos) && (s5.find_first_not_of("0123456789") == std::string::npos)) {
        opt::samH=false;
        samSeq=hseq;
        return 2;
    }
    return 3;
}

inline void ntComp(const uint64_t hVal, uint16_t *t_Counter) {
    uint64_t indBit=opt::nSamp;
    if(hVal>>(63-opt::sBits) == 1) indBit=0;
    if(hVal>>(64-opt::sBits) == opt::sMask) indBit=1;
    if(indBit < opt::nSamp) {
        size_t shVal=hVal&(opt::rBuck-1);
        #pragma omp atomic
        ++t_Counter[indBit*opt::rBuck+shVal];
    }

}

inline void ntRead(const string &seq, uint16_t *t_Counter, size_t &totKmer) {
    ntHashIterator itr(seq, 1, opt::kmLen);
    while (itr != itr.end()) {
        ntComp((*itr)[0],t_Counter);
        ++itr;
        ++totKmer;
    }
}

void getEfq(std::ifstream &in, uint16_t *t_Counter, size_t &totKmer) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = static_cast<bool>(getline(in, seq));
        good = static_cast<bool>(getline(in, hseq));
        good = static_cast<bool>(getline(in, hseq));
        if(good)
            ntRead(seq, t_Counter, totKmer);
        good = static_cast<bool>(getline(in, hseq));
    }
}

void getEfa(std::ifstream &in, uint16_t *t_Counter, size_t &totKmer) {
    bool good = true;
    for(string seq, hseq; good;) {
        string line;
        good = static_cast<bool>(getline(in, seq));
        while(good&&seq[0]!='>') {
            line+=seq;
            good = static_cast<bool>(getline(in, seq));
        }
        ntRead(line, t_Counter, totKmer);
    }
}

void getEsm(std::ifstream &in, const std::string &samSeq, uint16_t *t_Counter, size_t &totKmer) {
    std::string samLine,seq;
    std::string s1,s2,s3,s4,s5,s6,s7,s8,s9,s11;
    if(opt::samH) {
        while(getline(in,samLine))
            if (samLine[0]!='@') break;
    }
    else
        samLine=samSeq;
    do {
        std::istringstream iss(samLine);
        iss>>s1>>s2>>s3>>s4>>s5>>s6>>s7>>s8>>s9>>seq>>s11;
        ntRead(seq, t_Counter, totKmer);
    } while(getline(in,samLine));
}

void compEst(const uint16_t *t_Counter, double &F0Mean, double fMean[]) {
    unsigned p[opt::nSamp][65536];
    for(size_t i=0; i<opt::nSamp; i++)
        for(size_t j=0; j<65536; j++)
            p[i][j]=0;

    for(size_t i=0; i<opt::nSamp; i++)
        for(size_t j=0; j<opt::rBuck; j++)
            ++p[i][t_Counter[i*opt::rBuck+j]];

    double pMean[65536];
    for(size_t i=0; i<65536; i++) pMean[i]=0.0;
    for(size_t i=0; i<65536; i++) {
        for(size_t j=0; j<opt::nSamp; j++)
            pMean[i]+=p[j][i];
        pMean[i] /= 1.0*opt::nSamp;
    }

    F0Mean= (ssize_t)((opt::rBits*log(2)-log(pMean[0])) * 1.0* ((size_t)1<<(opt::sBits+opt::rBits)));
    for(size_t i=0; i<65536; i++) fMean[i]=0;
    fMean[1]= -1.0*pMean[1]/(pMean[0]*(log(pMean[0])-opt::rBits*log(2)));
    for(size_t i=2; i<65536; i++) {
        double sum=0.0;
        for(size_t j=1; j<i; j++)
            sum+=j*pMean[i-j]*fMean[j];
        fMean[i]=-1.0*pMean[i]/(pMean[0]*(log(pMean[0])-opt::rBits*log(2)))-sum/(i*pMean[0]);
    }
    for(size_t i=1; i<65536; i++)
        fMean[i]=abs((ssize_t)(fMean[i]*F0Mean));
}

void getHist(const vector<string> &inFiles, const unsigned kLen, const unsigned nThr, size_t histArray[]) {

    double sTime = omp_get_wtime();

    opt::nThrd = nThr;
    opt::kmLen = kLen;
    
    size_t totalSize=0;
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i)
        totalSize += getInf(inFiles[file_i].c_str());
    if(totalSize<50000000000) opt::sBits=7;

    size_t totalKmers = 0;

    opt::rBuck = ((size_t)1) << opt::rBits;
    opt::sMask = (((size_t)1) << (opt::sBits-1))-1;
    uint16_t *t_Counter = new uint16_t [opt::nSamp*opt::rBuck]();

#ifdef _OPENMP
    omp_set_num_threads(opt::nThrd);
#endif

    #pragma omp parallel for schedule(dynamic)
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
        size_t totKmer=0;
        std::ifstream in(inFiles[file_i].c_str());
        std::string samSeq;
        unsigned ftype = getftype(in,samSeq);
        if(ftype==0)
            getEfq(in, t_Counter, totKmer);
        else if (ftype==1)
            getEfa(in, t_Counter, totKmer);
        else if (ftype==2)
            getEsm(in, samSeq, t_Counter, totKmer);
        else {
            std::cerr << "Error in reading file: " << inFiles[file_i] << std::endl;
            exit(EXIT_FAILURE);
        }
        in.close();
        #pragma omp atomic
        totalKmers+=totKmer;
    }

    double F0Mean=0.0;
    double fMean[65536];
    compEst(t_Counter, F0Mean, fMean);
    histArray[0] = totalKmers;
    histArray[1] = (size_t)F0Mean;
    for(size_t i=2; i<= opt::covMax+1; i++)
        histArray[i] = (size_t)fMean[i-1];
    delete [] t_Counter;
    std::cerr << "Runtime(sec): " <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
}
#endif /* NTCARD_H_ */

