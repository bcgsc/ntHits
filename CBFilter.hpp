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

    CBFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize, const char * fPath):
        m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize) {
        m_filter = new unsigned char [m_size];
        ifstream myFile(fPath, ios::in | ios::binary);
        myFile.seekg (0, ios::beg);
        myFile.read ((char *)m_filter, m_size);
        myFile.close();
    }

    CBFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize, unsigned repCap):
    m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize), m_reCap(repCap) {
        m_filter = new unsigned char [m_size];
        for(size_t i = 0; i < m_size; i++)
            m_filter[i]=0;
    }

    
    bool insert_and_test(const uint64_t *hVal) {
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            unsigned char old_byte = __sync_fetch_and_add(&m_filter[hLoc],1);
            if(old_byte <= m_reCap) return false;
        }
        return true;
    }

    void storeFilter(const char * fPath) const {
        ofstream myFile(fPath, ios::out | ios::binary);
        myFile.write(reinterpret_cast<char*>(m_filter), m_size);
        myFile.close();
    }

    unsigned getKmerSize() const {
        return m_kmerSize;
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
