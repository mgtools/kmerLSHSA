#include <chrono>
#include <iostream>
#include <map>
#include <mutex>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <thread>
#include <omp.h>

#include "../common/abundance.h"
#include "../hash/lshash.h"
#include "distance.h"
#include "funcAB.h"
#include "../io/ioMatrix.h"

//#include "../utils/threadpool.hpp"
//#include "../boost/bind.hpp"


//using namespace boost::threadpool;
using namespace std;
using namespace Core;
using namespace Hash;
using namespace Utility;

void p_lsh(vector<Abundance*> *hash_values, vector<int> *hash_keys, hashTable& hash_table, Abundance* abundance);

void merge_hashtable(vector<vector<Abundance*>>* final_table, int hash_func_num, const vector<vector<Abundance*>*>& part_hash_values,  const vector<vector<int>*>& part_hash_keys);

void free_hashtable(vector<vector<Abundance*>>* final_table, int s);

void free_abundance(vector<Abundance*> final_abundance);

void merge_abundance(vector<Abundance*>* final_abundance, const vector<vector<Abundance*>*>& part_abundance);

void p_cluster(vector<Abundance*> *part_ab, vector<Abundance*>* candidates, float threshold);

void nestedCluster(vector<Abundance*>* unknown_abundance_ptr, float similarity, int dim, int thread_to_use, bool verbose);

void Cluster(vector<Abundance*>* unknown_abundance_ptr, float min_similarity, float max_similarity, int cluster_iteration, unsigned int threads_to_use, int dim, int bucket_size_threshold, bool verbose);
