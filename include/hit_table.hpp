#ifndef HIT_TABLE_HPP
#define HIT_TABLE_HPP

#include <fstream>
#include <omp.h>
#include <stdint.h>
#include <string>

#include "utils.hpp"

namespace nthits {

class HitTable
{
private:
  struct TableEntry
  {
    std::string kmer;
    unsigned count;
  };

  const size_t table_size;
  TableEntry* entries;
  const unsigned lock_size;
  omp_lock_t* locks;

public:
  explicit HitTable(const size_t table_size, const unsigned lock_size = 65536)
    : table_size(table_size)
    , entries(new TableEntry[table_size])
    , lock_size(lock_size)
    , locks(new omp_lock_t[lock_size])
  {
    for (size_t i = 0; i < table_size; i++) {
      entries[i].count = 0;
    }
    for (unsigned i = 0; i < lock_size; i++) {
      omp_init_lock(&locks[i]);
    }
  }

  ~HitTable()
  {
    for (unsigned i = 0; i < lock_size; i++)
      omp_destroy_lock(&locks[i]);
    delete[] locks;
    delete[] entries;
  }

  void insert(const uint64_t hash_value, const std::string& kmer, unsigned n = 1)
  {
    std::string canonical = kmer;
    to_canonical(canonical);
    uint64_t i = 0, j;
    do {
      j = (hash_value + i) % table_size;
      if (entries[j].kmer == canonical) {
#pragma omp atomic
        entries[j].count += n;
      }
      ++i;
    } while (i != table_size && entries[j].count != 0);
    if (entries[j].count == 0) {
      omp_set_lock(&locks[(uint16_t)j]);
      entries[j].kmer = canonical;
      ++entries[j].count;
      omp_unset_lock(&locks[(uint16_t)j]);
    }
  }

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
  void remove(const uint64_t hash_value)
  { // TODO
  }
#pragma GCC diagnostic pop

  void save(const std::string& file_path, const unsigned min_count)
  {
    std::ofstream outFile(file_path);
    for (size_t i = 0; i < table_size; i++)
      if (entries[i].count != 0)
        outFile << entries[i].kmer << "\t" << entries[i].count + min_count << std::endl;
    outFile.close();
  }
};

}

#endif