#ifndef KLSH_KMC_READER_HPP
#define KLSH_KMC_READER_HPP

#include <iostream>
#include <vector>


#include "kmc_api/kmc_file.h"
#include "Kmer.h"
#include "../hash/HashTables.h"
#include "kmc_api/kmer_defs.h"
//#include "../utils/threadpool.hpp"
//#include <boost/bind.hpp>
//using namespace boost::threadpool;
using namespace std;

void KmcRead(string kmc_file_name, ckhmap_t *kmap_ptr, bool verbose,  unsigned int num_threads, size_t ksize);
//void KmcRead(string kmc_file_name, ckhmap_t *kmap_ptr, bool verbose, pool &tp, unsigned int num_threads);

float_t KmcCount(string kmc_file_name, ckhmap_t *kmap_ptr, bool verbose,  unsigned int num_threads, size_t ksize);
//float_t KmcCount(string kmc_file_name, ckhmap_t *kmap_ptr, bool verbose, pool &tp, unsigned int num_threads);

#endif //KLSH_KMC_READER_HPP
