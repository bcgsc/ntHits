#include "args.hpp"
#include "ntcard.hpp"
#include "utils.hpp"

#include <argparse/argparse.hpp>
#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <btllib/nthash.hpp>
#include <btllib/seq_reader.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using cbf_counter_t = uint8_t;

struct entry
{
	std::string kmer;
	unsigned count;
};

bool
hit_search_insert(
    const uint64_t hash_value,
    const std::string& kmer,
    omp_lock_t* locks,
    entry* table,
    size_t table_size)
{
	uint64_t i = 0, j;
	do {
		j = (hash_value + i) % table_size;
		if (table[j].kmer == kmer) {
#pragma omp atomic
			++table[j].count;
			return true;
		}
		++i;
	} while (i != table_size && table[j].count != 0);
	if (table[j].count == 0) {
		omp_set_lock(&locks[(uint16_t)j]);
		table[j].kmer = kmer;
		++table[j].count;
		omp_unset_lock(&locks[(uint16_t)j]);
		return false;
	}
	return false;
}

unsigned
hit_search(const uint64_t hash_value, const std::string& kmer, entry* table, size_t table_size)
{
	uint64_t i = 0, j;
	do {
		j = (hash_value + i) % table_size;
		if (table[j].kmer == kmer) {
			return table[j].count;
		}
		++i;
	} while (i != table_size && table[j].count != 0);
	return 0;
}

void
f_hit(
    const std::string& file_path,
    omp_lock_t* locks,
    btllib::KmerBloomFilter& distincts,
    btllib::CountingBloomFilter<cbf_counter_t>& cbf,
    entry* table,
    ProgramArguments& args)
{
	unsigned seq_reader_mode;
	if (args.long_mode()) {
		seq_reader_mode = btllib::SeqReader::Flag::LONG_MODE;
	} else {
		seq_reader_mode = btllib::SeqReader::Flag::SHORT_MODE;
	}
	btllib::SeqReader reader(file_path, seq_reader_mode);
#pragma omp parallel shared(reader)
	for (const auto& record : reader) {
		btllib::NtHash nth(record.seq, args.get_num_hashes() + 1, args.get_kmer_length());
		while (nth.roll()) {
			if (!distincts.contains_insert(nth.hashes())) {
				std::string canonKmer = record.seq.substr(nth.get_pos(), args.get_kmer_length());
				getCanon(canonKmer);
				if (args.get_hit_cap() > 1) {
					if (cbf.insert_thresh_contains(nth.hashes(), args.get_hit_cap() - 1)) {
						hit_search_insert(
						    nth.hashes()[0], canonKmer, locks, table, args.get_hit_size());
					}
				} else {
					hit_search_insert(
					    nth.hashes()[0], canonKmer, locks, table, args.get_hit_size());
				}
			}
		}
	}
}

void
b_hit(
    const std::string& file_path,
    btllib::KmerBloomFilter& mydBF,
    btllib::CountingBloomFilter<cbf_counter_t>& mycBF,
    btllib::KmerBloomFilter& myhBF,
    ProgramArguments& args)
{
	unsigned seq_reader_mode;
	if (ProgramArguments::get_instance().long_mode()) {
		seq_reader_mode = btllib::SeqReader::Flag::LONG_MODE;
	} else {
		seq_reader_mode = btllib::SeqReader::Flag::SHORT_MODE;
	}
	btllib::SeqReader reader(file_path, seq_reader_mode);
#pragma omp parallel shared(reader)
	for (const auto record : reader) {
		btllib::NtHash nth(record.seq, args.get_num_hashes() + 1, args.get_kmer_length());
		while (nth.roll()) {
			if (mydBF.contains_insert(nth.hashes())) {
				if (args.get_hit_cap() > 1) {
					if (mycBF.insert_thresh_contains(nth.hashes(), args.get_hit_cap() - 1))
						myhBF.insert(nth.hashes());
				} else {
					myhBF.insert(nth.hashes());
				}
			}
		}
	}
}

int
main(int argc, char** argv)
{
	double sTime = omp_get_wtime();

	auto args = ProgramArguments::get_instance();
	args.parse(argc, argv);

	if (args.use_ntcard()) {
		size_t histArray[10002];
		getHist(args.get_input_files(), args.get_kmer_length(), args.get_num_threads(), histArray);

		int histIndex = 2, errCov = 1;
		while (histIndex <= 10000 && histArray[histIndex] > histArray[histIndex + 1])
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

		if (args.solid()) {
			if (args.get_hit_cap() == 0)
				args.set_hit_cap(minCov);
			std::cerr << "Errors k-mer coverage: " << args.get_hit_cap() << std::endl;
		}

		if (!args.solid()) {
			if (args.get_hit_cap() == 0)
				args.set_hit_cap(1.75 * maxCov);
			std::cerr << "Errors k-mer coverage: " << minCov << std::endl;
			std::cerr << "Median k-mer coverage: " << maxCov << std::endl;
			std::cerr << "Repeat k-mer coverage: " << args.get_hit_cap() << std::endl;
		}

		args.set_dbf_size(args.get_bits() * histArray[1]);
		args.set_cbf_size(args.get_bytes() * (histArray[1] - histArray[2]));
		size_t hitCount = histArray[1];
		for (unsigned i = 2; i <= args.get_hit_cap() + 1; i++)
			hitCount -= histArray[i];
		args.set_hit_size(args.out_bloom() ? hitCount * args.get_m() : hitCount * 3);

		std::cerr << "Approximate# of distinct k-mers: " << histArray[1] << "\n";
		std::cerr << "Approximate# of solid k-mers: " << hitCount << "\n";
	} else {
		args.set_dbf_size(args.get_bits() * args.get_f0());
		args.set_cbf_size(args.get_bytes() * (args.get_f0() - args.get_f1()));
		args.set_hit_size(args.get_m() * args.get_fr());
		std::cerr << "Approximate# of distinct k-mers: " << args.get_f0() << "\n";
		std::cerr << "Approximate# of solid k-mers: " << args.get_fr() << "\n";
	}

#ifdef _OPENMP
	omp_set_num_threads(args.get_num_threads());
#endif

	btllib::KmerBloomFilter mydBF(args.get_dbf_size() / 8, 3, args.get_kmer_length());
	btllib::CountingBloomFilter<cbf_counter_t> mycBF(args.get_cbf_size(), args.get_num_hashes());

	if (args.out_bloom()) {
		btllib::KmerBloomFilter myhBF(
		    args.get_hit_size() / 8, args.get_num_hashes() + 1, args.get_kmer_length());
		for (const auto& file_path : args.get_input_files()) {
			b_hit(file_path, mydBF, mycBF, myhBF, args);
		}
		std::string hbf_out_path;
		if (args.get_prefix().empty() && args.solid()) {
			hbf_out_path = "solids_k" + std::to_string(args.get_kmer_length()) + ".bf";
		} else if (args.get_prefix().empty() && !args.solid()) {
			hbf_out_path = "repeat_k" + std::to_string(args.get_kmer_length()) + ".bf";
		} else {
			hbf_out_path =
			    args.get_prefix() + "_k" + std::to_string(args.get_kmer_length()) + ".bf";
		}
		myhBF.save(hbf_out_path);
	} else {
		entry* hitTable = new entry[args.get_hit_size()];
		for (size_t i = 0; i < args.get_hit_size(); i++)
			hitTable[i].count = 0;
		const unsigned lockSize = 65536;
		omp_lock_t* locks = new omp_lock_t[lockSize];

		for (unsigned i = 0; i < lockSize; i++)
			omp_init_lock(&locks[i]);

		for (const auto& file_path : args.get_input_files()) {
			f_hit(file_path, locks, mydBF, mycBF, hitTable, args);
		}

		for (unsigned i = 0; i < lockSize; i++)
			omp_destroy_lock(&locks[i]);
		delete[] locks;

		std::stringstream hstm;
		if (args.get_prefix().empty()) {
			if (args.solid())
				hstm << "solids_k" << args.get_kmer_length() << ".rep";
			else
				hstm << "repeat_k" << args.get_kmer_length() << ".rep";
		} else
			hstm << args.get_prefix() << "_k" << args.get_kmer_length() << ".rep";
		std::ofstream outFile(hstm.str().c_str());
		for (size_t i = 0; i < args.get_hit_size(); i++)
			if (hitTable[i].count != 0)
				outFile << hitTable[i].kmer << "\t" << hitTable[i].count + args.get_hit_cap()
				        << "\n";
		outFile.close();

		delete[] hitTable;
	}

	std::cerr << "Total time for computing repeat content in (sec): " << std::setprecision(4)
	          << std::fixed << omp_get_wtime() - sTime << "\n";
	return 0;
}
