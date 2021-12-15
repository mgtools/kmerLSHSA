#ifndef __UTILITY_IO_MAT__
#define __UTILITY_IO_MAT__

#include <algorithm>
#include <cassert>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <chrono>
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
#include <unistd.h>
#include <utility>
#include <vector>

#include "../common/abundance.h"
#include "../function/funcAB.h"

using namespace std;
using namespace Core;
using namespace std::chrono;
namespace Utility{
class IOMat {
 public:
   static int parseLine(char* line);

   static int getValue(int i);

   static float convert(char const* source, char ** endPtr ) ;

   static void ReadMatrix(vector<Abundance*>* AbundanceMat, bool normalization,  string file_name ) ;

   static void ReadClusterAll(vector<Abundance*>* AbundanceMat, int num_samples,  string file_name, bool verbose );

   static void ReadCluster(vector<Abundance*>* AbundanceMat, int num_samples,  string file_name, streamoff start_line, uint64_t num_lines, bool verbose);

   static void SaveResult(const vector<Abundance*>* abs_all, string out_dir, bool delfile, int ignore_small, bool verbose);

   static void SaveMatrix(const vector<Abundance*>* abs_all, string out_file, bool delfile, int ignore_small);

   static void SaveBinary(const vector<Abundance*>* abs_all, string out_file, bool delfile, int ignore_small, bool verbose);

   static void convertHTMat(uint16_t **ary_count, vector<float_t> &v_kmers, int tot_sample, bool verbose, uint64_t batch_size, streamoff batch_offset,  vector<Abundance*>* unknown_abundance_ptr);
};

}  // namespace Utility
#endif
