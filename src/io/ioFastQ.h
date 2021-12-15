#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>

#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <time.h>
#include <limits.h>
//#include <atomic>
#include <chrono>

#include <stdint.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <functional>

#include <getopt.h>

#include "../common/Common.h"
//#include "utils.hpp"
#include "../utils/fastq.h"
#include "../kmer/kmc_api/kmc_file.h"
#include "../kmer/kmc_reader.h"
#include "../kmer/Kmer.h"

//#include "../utils/threadpool.hpp"
//#include <boost/bind.hpp>
//#include <boost/thread.hpp>
using namespace std;
//using namespace boost::threadpool;
using namespace std::chrono;

namespace Utility{
class IOFQ{
    public:

//identification of distinctive reads using multi-threads
    static void CheckRead(uset_t *g_kmer_ptr, vector<ReadEntry> &read_vec, vector<int> &record_vec, unsigned int num_threads, int tid, float kmer_vote);

//extract distinctive reads for each sample
    static void ReadExtract(uset_t *g_kmer_ptr, vector<string> &files, string output, float kmer_vote, bool verbose, unsigned int num_threads);
    //static void ReadExtract(uset_t *g_kmer_ptr, vector<string> &files, string output, float kmer_vote, bool verbose, pool &tp, unsigned int num_threads);

    static void Extracting(vector<string> samples, uset_t *g_kmer_ptr, string out, int num_threads, float kmer_vote, bool verbose);
};
}

