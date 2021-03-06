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
#include "vendor/ntHash/nthash.hpp"

using namespace std;


class CBFilter {
public:

    CBFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize, unsigned repCap):
    m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize), m_reCap(repCap) {
        m_filter = new unsigned char [m_size]();
    }

    bool insert_and_test(const uint64_t *hVal) {
        bool greaterFlag = true;
        unsigned minCount = m_filter[hVal[0] % m_size];
        if(minCount < m_reCap)
            greaterFlag = false;
        for (unsigned i = 1; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            if(m_filter[hLoc] < m_reCap) {
                if(m_filter[hLoc] < minCount)
                    minCount = m_filter[hLoc];
                greaterFlag = false;
            }
        }
        if(!greaterFlag) {
            //unsigned curVal, newVal;
            for (unsigned i = 0; i < m_hashNum; i++) {
                size_t hLoc = hVal[i] % m_size;
                if(m_filter[hLoc] == minCount) {
                    //do {
                    //    curVal = m_filter[hLoc];
                    //    newVal = curVal + 1;
                    //    if (newVal < curVal)
                    //        break;
                    //} while(!__sync_bool_compare_and_swap(&m_filter[hLoc], curVal, newVal));
                    #pragma omp atomic
                    ++m_filter[hLoc];
                }
            }
        }
        return greaterFlag;
    }

    ~CBFilter() {
        delete[] m_filter;
    }

private:
    CBFilter(const CBFilter& that); //to prevent copy construction
    unsigned char *m_filter;
    size_t m_size;
    unsigned m_hashNum;
    unsigned m_kmerSize;
    unsigned m_reCap;
};

#endif /* CBFILTER_H_ */
