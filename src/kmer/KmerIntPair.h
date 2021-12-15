#ifndef KLSH_KMER_INT_PAIR_HPP
#define KLSH_KMER_INT_PAIR_HPP

#include "Kmer.h"

// TODO: switch out for a templated version
struct KmerIntPair {
  KmerIntPair() {};
  KmerIntPair(const Kmer &km, uint32_t k);

  //char v[sizeof(Kmer)+sizeof(uint32_t)];//sizeof(char) to sizeof(int)
  char v[sizeof(Kmer)];
  uint32_t cnt;
  uint32_t GetVal() const;
  void SetVal(const uint32_t k);
  const Kmer& GetKey() const;
  void SetKey(const Kmer& km);

  static const size_t KmerOffset = 0;
  static const size_t IntOffset = sizeof(Kmer);
};



struct SelectKmerKey {
  const Kmer& operator()(const KmerIntPair &p) const {
    return p.GetKey();
  }
};

struct SetKmerKey {
  void operator()(KmerIntPair *value, const Kmer& km);
};


#endif // SA_KMER_INT_PAIR_HPP
