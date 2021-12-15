#include <cstdlib>

#include "Kmer.h"
#include "KmerIntPair.h"


KmerIntPair::KmerIntPair(const Kmer &km, uint32_t val) {
  SetKey(km);
  SetVal(val);
}

void KmerIntPair::SetVal(const uint32_t val) {
  //char val8 = (val > 0xFF) ?  0xFF : (char)val;
  //memcpy(&this->v + KmerIntPair::IntOffset, &val8, sizeof(uint8_t));
  //this->v[KmerIntPair::IntOffset] = val8;
  //this->v[KmerIntPair::IntOffset] = val;
  KmerIntPair::cnt = val;
  //memcpy(this+KmerIntPair::IntOffset, &val, sizeof(uint32_t));

}

uint32_t KmerIntPair::GetVal() const {
  //uint8_t tmp = *reinterpret_cast<const uint8_t*>(this+KmerIntPair::IntOffset);
  //2/12/2015, uint8_t to uint32_t
  //uint32_t tmp = *reinterpret_cast<const uint32_t*>(this+KmerIntPair::IntOffset);
  //return tmp;
  //return (uint8_t)this->v[KmerIntPair::IntOffset];
  return KmerIntPair::cnt;
}

const Kmer& KmerIntPair::GetKey() const {
  return *reinterpret_cast<const Kmer*>(this + KmerIntPair::KmerOffset);
}

void KmerIntPair::SetKey(const Kmer& km) {
  memcpy(this, &km, sizeof(Kmer));
}

void SetKmerKey::operator()(KmerIntPair *value, const Kmer& km) {
  memcpy(value + KmerIntPair::KmerOffset, &km, sizeof(Kmer));
}


