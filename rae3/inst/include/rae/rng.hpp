// RAE 3.0 core -- self-contained random numbers.
// Own generator and samplers so that the R package and the CLI produce bit-identical results for a
// seed, independent of the C++ standard library in use.  [INFERRED G-11f]
#pragma once
#include <cmath>
#include <cstdint>

namespace rae {

class Rng {
  uint64_t s[4];
  static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  static uint64_t splitmix(uint64_t& z) {
    z += 0x9e3779b97f4a7c15ULL;
    uint64_t r = z;
    r = (r ^ (r >> 30)) * 0xbf58476d1ce4e5b9ULL;
    r = (r ^ (r >> 27)) * 0x94d049bb133111ebULL;
    return r ^ (r >> 31);
  }
  bool have_spare = false;
  double spare = 0;

 public:
  // stream: independent sub-stream (e.g. sample index, permutation round)
  explicit Rng(uint64_t seed, uint64_t stream = 0) {
    uint64_t z = seed * 0x9e3779b97f4a7c15ULL + stream * 0xd1b54a32d192ed03ULL + 0x2545f4914f6cdd1dULL;
    for (int i = 0; i < 4; ++i) s[i] = splitmix(z);
  }
  uint64_t next() {  // xoshiro256**
    const uint64_t result = rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl(s[3], 45);
    return result;
  }
  double unif() { return (double)(next() >> 11) * (1.0 / 9007199254740992.0); }  // [0, 1)
  uint32_t below(uint32_t n) {  // unbiased integer in [0, n)
    uint64_t m = (uint64_t)(uint32_t)(next() >> 32) * n;
    uint32_t l = (uint32_t)m;
    if (l < n) {
      uint32_t t = (0u - n) % n;
      while (l < t) { m = (uint64_t)(uint32_t)(next() >> 32) * n; l = (uint32_t)m; }
    }
    return (uint32_t)(m >> 32);
  }
  double norm() {  // Box-Muller (polar form)
    if (have_spare) { have_spare = false; return spare; }
    double u, v, q;
    do { u = 2 * unif() - 1; v = 2 * unif() - 1; q = u * u + v * v; } while (q >= 1 || q == 0);
    double f = std::sqrt(-2 * std::log(q) / q);
    spare = v * f; have_spare = true;
    return u * f;
  }
};

}  // namespace rae
