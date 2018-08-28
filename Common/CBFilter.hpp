/*
 *
 * CBFilter.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */


#ifndef CBFILTER_H_
#define CBFILTER_H_
#include <string>
#include <stdint.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <cstring>
#include "nthash.hpp"

using namespace std;


class CBFilter {
public:

    CBFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize, unsigned repCap):
    m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize), m_reCap(repCap) {
        m_filter = new unsigned char [m_size];
        for(size_t i = 0; i < m_size; i++) m_filter[i]=0;
        for(unsigned i = 0; i < 65536; i++) omp_init_lock(&locks1[i]);
        for(unsigned i = 0; i < 65536; i++) omp_init_lock(&locks2[i]);
    }

    bool insert_and_test(const uint64_t *hVal) {
        omp_set_lock(&locks1[(uint16_t)(hVal[0]%m_size)]);
        omp_set_lock(&locks2[(uint16_t)(hVal[1]%m_size)]);
        bool greaterFlag = true;
        unsigned minCount = 256;
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            if(m_filter[hLoc] < m_reCap) {
                if(m_filter[hLoc] < minCount)
                    minCount = m_filter[hLoc];
                greaterFlag = false;
            }
        }
        if(!greaterFlag) {
            for (unsigned i = 0; i < m_hashNum; i++) {
                size_t hLoc = hVal[i] % m_size;
                if(m_filter[hLoc] == minCount) {
                    //#pragma omp atomic
                    ++m_filter[hLoc];
                }
            }
        }
        omp_unset_lock(&locks1[(uint16_t)(hVal[0]%m_size)]);
        omp_unset_lock(&locks2[(uint16_t)(hVal[1]%m_size)]);
        return greaterFlag;
    }

    ~CBFilter() {
        delete[] m_filter;
        for(unsigned i = 0; i < 65536; i++) omp_destroy_lock(&locks1[i]);
        for(unsigned i = 0; i < 65536; i++) omp_destroy_lock(&locks2[i]);
    }

private:
    CBFilter(const CBFilter& that); //to prevent copy construction
    unsigned char *m_filter;
    omp_lock_t locks1[65536];
    omp_lock_t locks2[65536];
    size_t m_size;
    unsigned m_hashNum;
    unsigned m_kmerSize;
    unsigned m_reCap;
};

#endif /* CBFILTER_H_ */
