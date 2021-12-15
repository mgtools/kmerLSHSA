#ifndef SA_HASHTABLES_HPP
#define SA_HASHTABLES_HPP

#include<unordered_set>
//#include "google/sparse_hash_map"
//#include "google/sparsehash/sparsehashtable.h"

#include "../utils/libcuckoo/cuckoohash_map.hh"

#include "hash.h"
#include "../kmer/KmerIntPair.h"


//using google::sparse_hash_map;

//typedef google::sparse_hashtable<KmerIntPair, Kmer, KmerHash, SelectKmerKey, SetKmerKey, std::equal_to<Kmer>, std::allocator<KmerIntPair> > hmap_t;

//typedef google::sparse_hash_map<Kmer, float, KmerHash> hmapq_t;

typedef std::unordered_set<Kmer, KmerHash, std::equal_to<Kmer>, std::allocator<Kmer> > uset_t;

typedef cuckoohash_map<Kmer, uint16_t, KmerHash, std::equal_to<Kmer> > ckhmap_t;
//


#endif // SA_HASHTABLES_HPP
