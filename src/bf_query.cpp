#define PROGRAM_NAME "ntHits-BFQ"
#define PROGRAM_VERSION nthits::VERSION
#define PROGRAM_DESCRIPTION "Query tool for ntHits' output Bloom filter"
#define PROGRAM_COPYRIGHT "Copyright 2022 Canada's Michael Smith Genome Science Centre"

#include <argparse/argparse.hpp>
#include <btllib/bloom_filter.hpp>
#include <btllib/counting_bloom_filter.hpp>
#include <btllib/seq_reader.hpp>
#include <iostream>
#include <vector>

#include "nthits.hpp"

struct ProgramArguments
{
  std::string bf_path, ex_path, data_path;
  bool is_cbf;
};

ProgramArguments
parse_args(int argc, char** argv)
{
  ProgramArguments args;

  argparse::ArgumentParser parser(PROGRAM_NAME, PROGRAM_VERSION);
  parser.add_description(PROGRAM_DESCRIPTION);
  parser.add_epilog(PROGRAM_COPYRIGHT);

  parser.add_argument("--cbf")
    .help("Treat input file as a counting Bloom filter and output k-mer counts")
    .default_value(false)
    .implicit_value(true);

  parser.add_argument("-ex").help("Path to the Bloom filter containing excluded k-mers");

  parser.add_argument("bf_path").help("Input Bloom filter file").required();
  parser.add_argument("data_path")
    .help("Path to input dataset, use '-' to read from stdin.")
    .required();

  try {
    parser.parse_args(argc, argv);
  } catch (const std::runtime_error& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << parser;
    std::exit(1);
  }

  args.bf_path = parser.get("bf_path");
  args.data_path = parser.get("data_path");
  args.is_cbf = parser.get<bool>("--cbf");
  if (parser.is_used("-ex")) {
    args.ex_path = parser.get("-ex");
  }

  return args;
}

#define QUERY_CODE(BF, OUT)                                                                        \
  std::cerr << "Loading Bloom filter... " << std::flush;                                           \
  BF std::cerr << "DONE" << std::endl;                                                             \
  btllib::SeqReader reader(args.data_path, btllib::SeqReader::Flag::SHORT_MODE);                   \
  for (const auto& record : reader) {                                                              \
    btllib::NtHash h(record.seq, bf.get_hash_num(), bf.get_k());                                   \
    while (h.roll()) {                                                                             \
      kmer = record.seq.substr(h.get_pos(), bf.get_k());                                           \
      std::cout << kmer << "\t\t" << (OUT) << std::endl;                                           \
    }                                                                                              \
  }

int
main(int argc, char** argv)
{
  auto args = parse_args(argc, argv);

  std::string kmer, seq;

  if (args.is_cbf) {
    QUERY_CODE(btllib::KmerCountingBloomFilter8 bf(args.bf_path);, bf.contains(h.hashes()))
  } else if (args.ex_path.empty()) {
    QUERY_CODE(btllib::KmerBloomFilter bf(args.bf_path);
               , bf.contains(h.hashes()) > 0 ? "YES" : "NO")
  } else {
    QUERY_CODE(btllib::KmerBloomFilter bf(args.bf_path);
               btllib::KmerBloomFilter bf_ex(args.ex_path);
               , bf.contains(h.hashes()) > 0 && bf_ex.contains(h.hashes()) == 0 ? "YES" : "NO")
  }

  return 0;
}