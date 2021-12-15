#ifndef __CORE_AB_H__
#define __CORE_AB_H__

#include <algorithm>
#include <cassert>
#include <cstring>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <map>
#include <random>
#include <queue>
#include <set>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <unistd.h>
#include <utility>
#include <vector>

#include "../common/abundance.h"

#include "../utils/alglib-3.17.0/src/statistics.h"
#include "../utils/alglib-3.17.0/src/ap.h"

using namespace std;
using namespace Core;

namespace Core {
class AB {
 public:
  static void randAbundance(Abundance* abundance, const Abundance& ab, float scale);
  static void SetAbundance(Abundance* abundance,const vector<uint64_t>& ids, const vector<float>& values) ;
  static void SetNullAbundance(Abundance* abundance);
  static void SetConsensus(Abundance* abundance, const Abundance& ab1, const Abundance& ab2);
  static void WRS(unordered_set<uint64_t> *group1, unordered_set<uint64_t> *group2, Abundance* abundance, int num_sample1, int num_sample2, float pvalue_thresh, int size_thresh);
  static void dynamicST(unordered_set<uint64_t> *group1, unordered_set<uint64_t> *group2, vector<Abundance*>* clusteredab_ptr, int num_sample1, int num_sample2, float pvalue_thresh, int size_thresh, uint64_t extracted_size_thresh, bool verbose);
  static void SetMean(Abundance* abundance,const vector<Abundance*>& local_ab);
  static bool isSameAb(Abundance* ab1, Abundance* ab2);
  static float convert(char const* source, char ** endPtr);
};

}  // namespace Utility
#endif
