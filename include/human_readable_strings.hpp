#ifndef HUMAN_READABLE_STRINGS_HPP
#define HUMAN_READABLE_STRINGS_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>

using Uint = uint64_t;

// taken from
// https://stackoverflow.com/questions/63411054/how-can-you-quickly-compute-the-integer-logarithm-for-any-base/63411055#63411055

constexpr unsigned
log2floor(uint64_t x)
{
  // implementation for C++17 using clang or gcc
  return x ? 63 - __builtin_clzll(x) : 0;
}

template<typename Uint>
constexpr Uint
logFloor_naive(Uint val, unsigned base)
{
  Uint result = 0;
  while (val /= base) {
    ++result;
  }
  return result;
}

template<typename Uint, size_t BASE>
constexpr std::array<uint8_t, std::numeric_limits<Uint>::digits>
makeGuessTable()
{
  decltype(makeGuessTable<Uint, BASE>()) result{};
  for (size_t i = 0; i < result.size(); ++i) {
    Uint pow2 = static_cast<Uint>(Uint{ 1 } << i);
    result.data()[i] = logFloor_naive(pow2, BASE);
  }
  return result;
}

// The maximum possible exponent for a given base that can still be represented
// by a given integer type.
// Example: maxExp<uint8_t, 10> = 2, because 10^2 is representable by an 8-bit unsigned
// integer but 10^3 isn't.
template<typename Uint, unsigned BASE>
constexpr Uint maxExp = logFloor_naive<Uint>(static_cast<Uint>(~Uint{ 0u }), BASE);

// the size of the table is maxPow<Uint, BASE> + 2 because we need to store the maximum power
// +1 because we need to contain it, we are dealing with a size, not an index
// +1 again because for narrow integers, we access guess+1
template<typename Uint, size_t BASE>
constexpr std::array<uint64_t, maxExp<Uint, BASE> + 2>
makePowerTable()
{
  decltype(makePowerTable<Uint, BASE>()) result{};
  uint64_t x = 1;
  for (size_t i = 0; i < result.size(); ++i, x *= BASE) {
    result.data()[i] = x;
  }
  return result;
}

// If our base is a power of 2, we can convert between the
// logarithms of different bases without losing any precision.
constexpr bool
isPow2or0(uint64_t val)
{
  return (val & (val - 1)) == 0;
}

template<size_t BASE = 10, typename Uint>
constexpr Uint
logFloor(Uint val)
{
  if constexpr (isPow2or0(BASE)) {
    return log2floor(val) / log2floor(BASE);
  } else {
    constexpr auto guesses = makeGuessTable<Uint, BASE>();
    constexpr auto powers = makePowerTable<Uint, BASE>();

    uint8_t guess = guesses[log2floor(val)];

    // Accessing guess + 1 isn't always safe for 64-bit integers.
    // This is why we need this condition. See below for more details.
    if constexpr (sizeof(Uint) < sizeof(uint64_t) || guesses.back() + 2 < (int)powers.size()) {
      return guess + (val >= powers[guess + 1]);
    } else {
      return guess + (val / BASE >= powers[guess]);
    }
  }
}

// taken from
// https://stackoverflow.com/questions/63511627/how-can-i-stringify-a-fraction-with-n-decimals-in-c/63511628#63511628

std::string
stringifyFraction(const uint64_t num, const unsigned den, const unsigned precision)
{
  constexpr unsigned base = 10;

  // prevent division by zero if necessary
  if (den == 0) {
    return "inf";
  }

  // integral part can be computed using regular division
  std::string result = std::to_string(num / den);

  // perform first step of long division
  // also cancel early if there is no fractional part
  unsigned tmp = num % den;
  if (tmp == 0 || precision == 0) {
    return result;
  }

  // reserve characters to avoid unnecessary re-allocation
  result.reserve(result.size() + precision + 1);

  // fractional part can be computed using long division
  result += '.';
  for (size_t i = 0; i < precision; ++i) {
    tmp *= base;
    char nextDigit = '0' + static_cast<char>(tmp / den);
    result.push_back(nextDigit);
    tmp %= den;
  }

  return result;
}

// taken from
// https://stackoverflow.com/questions/63512258/how-can-i-print-a-human-readable-file-size-in-c-without-a-loop

// use SFINAE to only allow base 1000 or 1024
template<size_t BASE = 1000, std::enable_if_t<BASE == 1000 || BASE == 1024, int> = 0>
std::string
stringifyFileSize(uint64_t size, unsigned precision = 1) noexcept
{
  constexpr const char FILE_SIZE_UNITS[8][3]{ "B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB" };

  // The linked post about computing the integer logarithm
  // explains how to compute this.
  // This is equivalent to making a table: {1, 1000, 1000 * 1000, ...}
  // or {1, 1024, 1024 * 1024, ...}
  constexpr auto powers = makePowerTable<Uint, BASE>();

  unsigned unit = logFloor<BASE>(size);

  // Your numerator is size, your denominator is 1000^unit or 1024^unit.
  std::string result = stringifyFraction(size, powers[unit], precision);
  result.reserve(result.size() + 5);

  // Optional: Space separating number from unit. (usually looks better)
  result.push_back(' ');
  char first = FILE_SIZE_UNITS[unit][0];
  result.push_back(first);

  // Don't insert anything more in case of single bytes.
  if (unit != 0) {
    if constexpr (BASE == 1024) {
      result.push_back('i');
    }
    result.push_back(FILE_SIZE_UNITS[unit][1]);
  }

  return result;
}

#endif