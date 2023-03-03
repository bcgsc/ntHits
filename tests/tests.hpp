#ifndef NTHITS_TESTS_HPP
#define NTHITS_TESTS_HPP

#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <btllib/seq_reader.hpp>
#include <btllib/util.hpp>
#include <string>
#include <unordered_map>
#include <vector>

#include "hit_table.hpp"

#define BF_SIZE 1024 * 1024
#define NUM_HASHES 3

const std::string data_file = btllib::get_dirname(__FILE__) + "/data.fa";
const std::string freq_file = btllib::get_dirname(__FILE__) + "/freq.hist";
const std::string counts_file = btllib::get_dirname(__FILE__) + "/counts.txt";

std::vector<std::string>
load_test_data()
{
  std::vector<std::string> data;
  btllib::SeqReader reader(data_file, btllib::SeqReader::Flag::SHORT_MODE);
  for (const auto& record : reader) {
    data.push_back(record.seq);
  }
  return data;
}

class Validator
{

private:
  std::vector<std::string>& data;
  std::unordered_map<std::string, unsigned> true_counts;
  unsigned min_count, max_count;

  void load_counts()
  {
    std::ifstream cfile(counts_file);
    std::string kmer;
    unsigned count;
    while (cfile >> kmer >> count) {
      true_counts[kmer] = count;
    }
  }

public:
  Validator(std::vector<std::string>& data, unsigned min_count = 0, unsigned max_count = 0)
    : data(data)
    , min_count(min_count)
    , max_count(max_count)
  {
    load_counts();
  }

  unsigned get_kmer_length() const { return true_counts.begin()->first.size(); }

  bool check_kmer_bf(btllib::KmerBloomFilter& bf)
  {
    for (auto& it : true_counts) {
      bool contains = bf.contains(it.first) == 1;
      bool should_contain = it.second >= min_count && (max_count == 0 || it.second <= max_count);
      if (should_contain != contains) {
        std::cout << it.first << " " << it.second << " " << contains << " " << should_contain
                  << std::endl;
        return false;
      }
    }
    return true;
  }

  bool check_kmer_cbf(btllib::KmerCountingBloomFilter8& cbf)
  {
    for (auto& it : true_counts) {
      bool contains = cbf.contains(it.first) > 0;
      bool should_contain = it.second >= min_count && (max_count == 0 || it.second <= max_count);
      if (should_contain != contains) {
        return false;
      }
    }
    return true;
  }

};

#endif