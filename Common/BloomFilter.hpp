/*
 *
 * BloomFilter.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */


#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_

#include <math.h>

#include "iostream"
#include "fstream"
#include <cstring>
#include "nthash.hpp"
#include "vendor/cpptoml/include/cpptoml.h"

using namespace std;

inline unsigned popCnt(unsigned char x) {
    return ((0x876543210 >>
             (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >>
            ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2))
    & 0xf;
}

class BloomFilter {
public:

	BloomFilter(const string& filterFilePath)
	  : m_filter(NULL)
	{
		loadFilter(filterFilePath);
	}

	void loadFilter(const string& filterFilePath)
	{
		std::ifstream file(filterFilePath);
		loadHeader(file);
		// NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
		file.read(reinterpret_cast<char*>(m_filter), (m_size + 7)/8);
		file.close();
	}

	void loadHeader(std::istream& file)
	{
		std::string magic_header(MAGIC_HEADER_STRING);
		(magic_header.insert(0, "[")).append("]");
		std::string line;
		std::getline(file, line);
		if (line != magic_header) {
			std::cerr << "ERROR: magic string does not match (likely version mismatch)\n"
			          << "Your magic string:                " << line << "\n"
			          << "CountingBloomFilter magic string: " << magic_header << std::endl;
			exit(EXIT_FAILURE);
		}
    }

    BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize):
        m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize) {
        m_filter = new unsigned char [(m_size + 7)/8]();
    }

    bool insert_make_change(const uint64_t *hVal) {
        bool change=false;
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            unsigned char old_byte = __sync_fetch_and_or(&m_filter[hLoc/8],(1<<(7-hLoc%8)));
            if((old_byte&(1<<(7-hLoc%8)))==0) change=true;
        }
        return change;
    }

    void insert(const uint64_t *hVal) {
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            __sync_or_and_fetch(&m_filter[hLoc / 8], (1 << (7 - hLoc % 8)));
        }
    }

    void insert(const char* kmer) {
	uint64_t hVal = NTC64(kmer, m_kmerSize);
        for (unsigned i = 0; i < m_hashNum; i++) {
	    uint64_t mhVal=NTE64(hVal, m_kmerSize, i);
            size_t hLoc = mhVal % m_size;
            __sync_or_and_fetch(&m_filter[hLoc / 8], (1 << (7 - hLoc % 8)));
        }
    }

    bool contains(const char* kmer) const {
        uint64_t hVal = NTC64(kmer, m_kmerSize);
        for (unsigned i = 0; i < m_hashNum; i++) {
            uint64_t mhVal=NTE64(hVal, m_kmerSize, i);
            size_t hLoc = mhVal % m_size;
            if ((m_filter[hLoc / 8] & (1 << (7 - hLoc % 8))) == 0)
                return false;
        }
        return true;
    }

	void writeHeader(std::ostream& out) const
	{
		/* Initialize cpptoml root table
		   Note: Tables and fields are unordered
		   Ordering of table is maintained by directing the table
		   to the output stream immediately after completion  */
		std::shared_ptr<cpptoml::table> root = cpptoml::make_table();

		/* Initialize bloom filter section and insert fields
		   and output to ostream */
		auto header = cpptoml::make_table();
		header->insert("KmerSize", m_kmerSize);
		header->insert("HashNum", m_hashNum);
		header->insert("BloomFilterSize", m_size);
		header->insert("BloomFilterSizeInBytes", (m_size + 7)/8);
		header->insert("dFPR",pow(getPop()*1.0/m_size, m_hashNum));
		header->insert("nEntry", 0);
		header->insert("Entry", 0);
		std::string magic(MAGIC_HEADER_STRING);
		root->insert(magic, header);
		out << *root;

		// Output [HeaderEnd]\n to ostream to mark the end of the header
		out << "[HeaderEnd]\n";
}


	/** Serialize the Bloom filter to a stream */
	friend std::ostream& operator<<(std::ostream& out, const BloomFilter& bloom)
	{
		bloom.writeHeader(out);
		// NOLINTNEXTLINE(google-readability-casting)
		out.write(reinterpret_cast<char*>(bloom.m_filter), (bloom.m_size + 7)/8);
		return out;
	}


    void storeFilter(const char * fPath) const {
		std::ofstream ofs(fPath, std::ios::out | std::ios::binary);
		std::cerr << "Writing a " << (m_size + 7)/8 << " byte filter to " << fPath
		          << " on disk.\n";
		ofs << *this;
		ofs.flush();
        ofs.close();
    }

    size_t getPop() const {
        size_t i, popBF=0;
        #pragma omp parallel for reduction(+:popBF)
        for(i=0; i<(m_size + 7)/8; i++)
            popBF = popBF + popCnt(m_filter[i]);
        return popBF;
    }
	
    size_t getSize() const {
        return m_size;
    }
     
    ~BloomFilter() {
        delete[] m_filter;
    }

private:
    BloomFilter(const BloomFilter& that);
    unsigned char * m_filter;
    size_t m_size;
    unsigned m_hashNum;
    unsigned m_kmerSize;
    static constexpr const char* MAGIC_HEADER_STRING = "BTLBloomFilter_v1";
};

#endif /* BLOOMFILTER_H_ */
