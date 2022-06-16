#ifndef NTHITS_ARGS_HPP
#define NTHITS_ARGS_HPP

#define PROGRAM_NAME "ntHits"
#define PROGRAM_VERSION "0.0.1"
#define PROGRAM_DESCRIPTION "Reports the most frequent k-mers in input files."
#define PROGRAM_COPYRIGHT "Copyright 2019 Canada's Michael Smith Genome Science Centre"

#include <string>
#include <vector>

class ProgramArguments
{
  private:
	unsigned t;
	unsigned k;
	unsigned h;
	unsigned hit_cap = 0;
	unsigned bits = 7;
	unsigned bytes = 6;
	unsigned m = 16;
	unsigned dbf_size;
	unsigned cbf_size;
	unsigned hit_size;
	std::string prefix;
	unsigned F0, f1, fr;
	bool _out_bloom;
	bool _solid;
	bool _long_mode;
	bool _use_ntcard;
	std::vector<std::string> input_files;
	std::vector<std::string> seeds;

	ProgramArguments() = default;

  public:
	static ProgramArguments& get_instance()
	{
		static ProgramArguments instance;
		return instance;
	}

	void parse(int argc, char** argv);

	[[nodiscard]] unsigned get_num_threads() const { return t; }
	[[nodiscard]] unsigned get_kmer_length() const { return k; }
	[[nodiscard]] unsigned get_num_hashes() const { return h; }
	[[nodiscard]] unsigned get_hit_cap() const { return hit_cap; }
	[[nodiscard]] unsigned get_bits() const { return bits; }
	[[nodiscard]] unsigned get_bytes() const { return bytes; }
	[[nodiscard]] unsigned get_m() const { return m; }
	[[nodiscard]] unsigned get_dbf_size() const { return dbf_size; }
	[[nodiscard]] unsigned get_cbf_size() const { return cbf_size; }
	[[nodiscard]] unsigned get_hit_size() const { return hit_size; }
	[[nodiscard]] std::string get_prefix() const { return prefix; }
	[[nodiscard]] unsigned get_f0() const { return F0; }
	[[nodiscard]] unsigned get_f1() const { return f1; }
	[[nodiscard]] unsigned get_fr() const { return fr; }
	[[nodiscard]] bool out_bloom() const { return _out_bloom; }
	[[nodiscard]] bool solid() const { return _solid; }
	[[nodiscard]] bool long_mode() const { return _long_mode; }
	[[nodiscard]] bool use_ntcard() const { return _use_ntcard; }
	[[nodiscard]] std::vector<std::string> get_input_files() const { return input_files; }
	[[nodiscard]] std::vector<std::string> get_seeds() const { return seeds; }

	void set_hit_cap(const unsigned x) { hit_cap = x; }
	void set_dbf_size(const unsigned x) { dbf_size = x; }
	void set_cbf_size(const unsigned x) { cbf_size = x; }
	void set_hit_size(const unsigned x) { hit_size = x; }
};

#endif // NTHITS_ARGS_HPP
