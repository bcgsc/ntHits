#include "args.hpp"
#include "nthits.hpp"
#include "nthits_min.hpp"
#include "nthits_min_max.hpp"
#include "user_interface.hpp"
#include "utils.hpp"

#include <argparse/argparse.hpp>
#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <btllib/nthash.hpp>
#include <btllib/seq_reader.hpp>
#include <cmath>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

#define INIT(CONSTRUCTOR_CALLS)                                                                    \
  std::cout << "Initializing... " << std::flush;                                                   \
  timer.start();                                                                                   \
  CONSTRUCTOR_CALLS                                                                                \
  timer.stop();                                                                                    \
  timer.print_done();

#define PROCESS(NTHITS_CALL)                                                                       \
  std::cout << "Processing data... " << std::flush;                                                \
  timer.start();                                                                                   \
  for (const auto file : args.input_files) {                                                       \
    btllib::SeqReader reader(file, get_flag(args.long_mode));                                      \
    _Pragma("omp parallel shared(reader)") for (const auto& record : reader)                       \
    {                                                                                              \
      NTHITS_CALL                                                                                  \
    }                                                                                              \
  }                                                                                                \
  timer.stop();                                                                                    \
  timer.print_done();

#define SAVE(SAVE_CALL)                                                                            \
  std::cout << "Saving output... " << std::flush;                                                  \
  timer.start();                                                                                   \
  SAVE_CALL                                                                                        \
  timer.stop();                                                                                    \
  timer.print_done();

#define PRINT_BF_STATS(NAME, BF, MIN_VERBOSITY)                                                    \
  if (args.verbosity > MIN_VERBOSITY) {                                                            \
    std::cout << std::endl << NAME << " stats:" << std::endl;                                      \
    print_bloom_filter_stats(BF.get_fpr(), args.fpr, BF.get_occupancy());                          \
    std::cout << std::endl;                                                                        \
  }

inline unsigned
get_flag(bool long_mode)
{
  if (long_mode) {
    return btllib::SeqReader::Flag::LONG_MODE;
  } else {
    return btllib::SeqReader::Flag::SHORT_MODE;
  }
}

int
main(int argc, char** argv)
{
  print_logo();

  auto args = ProgramArguments(argc, argv);

  if (args.verbosity > 1) {
    args.print();
  }

  omp_set_num_threads(args.num_threads);

  auto hist = nthits::load_ntcard_histogram(args.histogram_path);
  size_t hit_count;
  unsigned given_hit_cap = args.min_count;
  nthits::get_thresholds(hist, args.solid, hit_count, args.min_count, args.max_count);
  bool hit_cap_changed = given_hit_cap != args.min_count;

  size_t bf_size, cbf_size, hit_size;
  if (args.min_count > 1) {
    bf_size = nthits::get_bf_size(hist[1], args.num_hashes, args.seeds.size(), sqrt(args.fpr)) / 8;
  } else {
    bf_size = 1;
  }
  if (args.min_count > 2 || args.has_max_count) {
    cbf_size =
      nthits::get_bf_size(hist[1] - hist[2], args.num_hashes, args.seeds.size(), sqrt(args.fpr));
  } else {
    cbf_size = 1;
  }
  size_t optimal_bytes =
    nthits::get_bf_size(hit_count, args.num_hashes, args.seeds.size(), args.fpr);
  if (args.out_type == OutputType::BLOOM_FILTER) {
    hit_size = optimal_bytes / 8;
  } else if (args.out_type == OutputType::COUNTING_BLOOM_FILTER) {
    hit_size = optimal_bytes;
  } else {
    hit_size = hit_count * 3;
  }
  if (args.verbosity > 0) {
    print_updated_params(hit_count,
                         hist[1],
                         bf_size,
                         args.min_count,
                         hit_cap_changed,
                         cbf_size,
                         args.out_is_filter(),
                         hit_size,
                         args.verbosity);
  }

  Timer timer;
  bool out_bf = args.out_type == OutputType::BLOOM_FILTER;
  bool out_cbf = args.out_type == OutputType::COUNTING_BLOOM_FILTER;
  bool out_table = args.out_type == OutputType::HIT_TABLE;
  bool min = args.has_min_count && !args.has_max_count;
  bool min_max = args.has_max_count;
  bool no_thresh = !args.has_min_count && !args.has_max_count;
  bool using_seeds = args.using_seeds();

  if (out_bf && no_thresh && !using_seeds) {
    INIT(btllib::KmerBloomFilter hits(hit_size, args.num_hashes, args.kmer_length);)
    PROCESS(nthits::find_hits(record.seq, args.kmer_length, hits);)
    PRINT_BF_STATS("Output BF", hits, 0)
    SAVE(hits.save(args.out_file);)
  }

  if (out_bf && min && !using_seeds) {
    INIT(btllib::BloomFilter bf(bf_size, args.num_hashes);
         btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
         btllib::KmerBloomFilter hits(hit_size, args.num_hashes, args.kmer_length);)
    PROCESS(nthits::find_hits(record.seq, args.kmer_length, args.min_count, bf, cbf, hits);)
    PRINT_BF_STATS("Distinct k-mers BF", bf, 1)
    PRINT_BF_STATS("Intermediate CBF", cbf, 1)
    PRINT_BF_STATS("Output BF", hits, 0)
    SAVE(hits.save(args.out_file);)
  }

  if (out_bf && !min && using_seeds) {
    INIT(btllib::BloomFilter hits(hit_size, args.num_hashes);)
    PROCESS({
      for (const auto& seed : args.seeds) {
        nthits::find_hits(record.seq, seed, hits);
      }
    })
    PRINT_BF_STATS("Output BF", hits, 0)
    SAVE(hits.save(args.out_file);)
  }

  if (out_bf && min && using_seeds) {
    INIT(btllib::BloomFilter bf(bf_size, args.num_hashes);
         btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
         btllib::BloomFilter hits(hit_size, args.num_hashes);)
    PROCESS({
      for (const auto& seed : args.seeds) {
        nthits::find_hits(record.seq, seed, args.min_count, bf, cbf, hits);
      }
    })
    PRINT_BF_STATS("Distinct seeds BF", bf, 1)
    PRINT_BF_STATS("Intermediate CBF", cbf, 1)
    PRINT_BF_STATS("Output BF", hits, 0)
    SAVE(hits.save(args.out_file);)
  }

  if (out_bf && min_max && !using_seeds) {
    INIT(btllib::BloomFilter bf(bf_size, args.num_hashes);
         btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
         btllib::KmerBloomFilter hits(hit_size, args.num_hashes, args.kmer_length);
         btllib::KmerBloomFilter excludes(hit_size, args.num_hashes, args.kmer_length);)
    PROCESS(
      nthits::find_hits(
        record.seq, args.kmer_length, args.min_count, args.max_count, bf, cbf, hits, excludes);)
    PRINT_BF_STATS("Distinct k-mers BF", bf, 1)
    PRINT_BF_STATS("Intermediate CBF", cbf, 1)
    PRINT_BF_STATS("Hit k-mers BF", hits, 0)
    PRINT_BF_STATS("Excluded k-mers BF", excludes, 0)
    SAVE(hits.save(args.out_file); excludes.save(args.out_file);)
  }

  if (out_bf && min_max && using_seeds) {
    INIT(btllib::BloomFilter bf(bf_size, args.num_hashes);
         btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
         btllib::BloomFilter hits(hit_size, args.num_hashes);
         btllib::BloomFilter excludes(hit_size, args.num_hashes);)
    PROCESS({
      for (const auto& seed : args.seeds) {
        nthits::find_hits(
          record.seq, seed, args.min_count, args.max_count, bf, cbf, hits, excludes);
      }
    })
    PRINT_BF_STATS("Distinct seeds BF", bf, 1)
    PRINT_BF_STATS("Intermediate CBF", cbf, 1)
    PRINT_BF_STATS("Hit seeds BF", hits, 0)
    PRINT_BF_STATS("Excluded seeds BF", excludes, 0)
    SAVE(hits.save(args.out_file); excludes.save(args.out_file);)
  }

  if (out_cbf && no_thresh && !using_seeds) {
    INIT(btllib::KmerCountingBloomFilter<nthits::cbf_counter_t> hits(
           hit_size, args.num_hashes, args.kmer_length);)
    PROCESS(nthits::find_hits(record.seq, args.kmer_length, hits);)
    PRINT_BF_STATS("Output BF", hits, 0)
    SAVE(hits.save(args.out_file);)
  }

  if (out_cbf && min && !using_seeds) {
    INIT(btllib::BloomFilter bf(bf_size, args.num_hashes);
         btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
         btllib::KmerCountingBloomFilter<nthits::cbf_counter_t> hits(
           hit_size, args.num_hashes, args.kmer_length);)
    PROCESS(nthits::find_hits(record.seq, args.kmer_length, args.min_count, bf, cbf, hits);)
    PRINT_BF_STATS("Distinct k-mers BF", bf, 1)
    PRINT_BF_STATS("Intermediate CBF", cbf, 1)
    PRINT_BF_STATS("Output BF", hits, 0)
    SAVE(hits.save(args.out_file);)
  }

  if (out_cbf && min_max && !using_seeds) {
    INIT(btllib::BloomFilter bf(bf_size, args.num_hashes);
         btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
         btllib::KmerCountingBloomFilter<nthits::cbf_counter_t> hits(
           hit_size, args.num_hashes, args.kmer_length);)
    PROCESS({
      nthits::find_hits(
        record.seq, args.kmer_length, args.min_count, args.max_count, bf, cbf, hits);
    })
    PRINT_BF_STATS("Distinct k-mers BF", bf, 1)
    PRINT_BF_STATS("Intermediate CBF", cbf, 1)
    PRINT_BF_STATS("Output BF", hits, 0)
    SAVE(hits.save(args.out_file);)
  }

  if (out_cbf && no_thresh && using_seeds) {
    INIT(btllib::CountingBloomFilter<nthits::cbf_counter_t> hits(hit_size, args.num_hashes);)
    PROCESS({
      for (const auto& seed : args.seeds) {
        nthits::find_hits(record.seq, seed, hits);
      }
    })
    PRINT_BF_STATS("Output BF", hits, 0)
    SAVE(hits.save(args.out_file);)
  }

  if (out_cbf && min && using_seeds) {
    INIT(btllib::BloomFilter bf(bf_size, args.num_hashes);
         btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
         btllib::CountingBloomFilter<nthits::cbf_counter_t> hits(hit_size, args.num_hashes);)
    PROCESS({
      for (const auto& seed : args.seeds) {
        nthits::find_hits(record.seq, seed, args.min_count, bf, cbf, hits);
      }
    })
    PRINT_BF_STATS("Distinct k-mers BF", bf, 1)
    PRINT_BF_STATS("Intermediate CBF", cbf, 1)
    PRINT_BF_STATS("Output BF", hits, 0)
    SAVE(hits.save(args.out_file);)
  }

  if (out_cbf && min_max && using_seeds) {
    INIT(btllib::BloomFilter bf(bf_size, args.num_hashes);
         btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
         btllib::CountingBloomFilter<nthits::cbf_counter_t> hits(hit_size, args.num_hashes);)
    PROCESS({
      for (const auto& seed : args.seeds) {
        nthits::find_hits(record.seq, seed, args.min_count, args.max_count, bf, cbf, hits);
      }
    })
    PRINT_BF_STATS("Distinct k-mers BF", bf, 1)
    PRINT_BF_STATS("Intermediate CBF", cbf, 1)
    PRINT_BF_STATS("Output BF", hits, 0)
    SAVE(hits.save(args.out_file);)
  }

  if (out_table && no_thresh) {
    INIT(nthits::HitTable hits(hit_size);)
    PROCESS(nthits::find_hits(record.seq, args.kmer_length, hits);)
    SAVE(hits.save(args.out_file, args.min_count);)
  }

  if (out_table && min) {
    INIT(btllib::BloomFilter bf(bf_size, args.num_hashes);
         btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
         nthits::HitTable hits(hit_size);)
    PROCESS(nthits::find_hits(record.seq, args.kmer_length, args.min_count, bf, cbf, hits);)
    PRINT_BF_STATS("Distinct k-mers BF", bf, 1)
    PRINT_BF_STATS("Intermediate CBF", cbf, 1)
    SAVE(hits.save(args.out_file, args.min_count);)
  }

  if (out_table && min_max) {
    INIT(btllib::BloomFilter bf(bf_size, args.num_hashes);
         btllib::CountingBloomFilter<nthits::cbf_counter_t> cbf(cbf_size, args.num_hashes);
         nthits::HitTable hits(hit_size);)
    PROCESS({
      nthits::find_hits(
        record.seq, args.kmer_length, args.min_count, args.max_count, bf, cbf, hits);
    })
    PRINT_BF_STATS("Distinct k-mers BF", bf, 1)
    PRINT_BF_STATS("Intermediate CBF", cbf, 1)
    SAVE(hits.save(args.out_file, args.min_count);)
  }

  return 0;
}
