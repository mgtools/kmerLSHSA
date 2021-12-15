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
#include <map>
#include <chrono>

#include <stdint.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <functional>

#include <getopt.h>

#include <mutex>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <thread>

#include "../common/Common.h"
#include "../common/abundance.h"

#include "../function/distance.h"
#include "../function/funcAB.h"
#include "../function/cluster.h"

#include "../io/ioMatrix.h"

#include "../hash/lshash.h"

//#include "../utils/threadpool.hpp"
//#include <boost/bind.hpp>

//#include <boost/thread.hpp>

using namespace std;
//using namespace boost;
//using namespace boost::threadpool;
using namespace std::chrono;
using namespace Core;
using namespace Utility;

struct HyperParams
{
  float min_similarity, max_similarity;
  float scale;
  int cluster_iteration;
  int max_memory;
  uint64_t batch_size, extracted_size;
  bool verbose;
  unsigned int threads_to_use;
  string input;
  string output;
};

void PrintHyperParams(const HyperParams& params) {
  cout << "************ kmers Cluster Params Setting ****************" << endl;
  cout << "cluster iteration: " << params.cluster_iteration << endl;
  cout << "group1 result prefix: " << params.output << endl;
  cout << "group1 input file: " << params.input << endl;
  cout << "min similarity: " << params.min_similarity << endl;
  cout << "max similarity: " << params.max_similarity << endl;
  cout << "threds to use: " << params.threads_to_use << endl;
  cout << "p-value threshold: " << params.pval_thresh << endl;
  cout << "********************************************************" << endl;
}

void kLSH_PrintUsage() {
  cerr << "kmers LSH "<< KLSH_VERSION << endl << endl;
  cerr << "Clustering of k-mers [from KMC library] " << endl << endl;
  cerr << "Usage: kmerLSH -i1 -i2 -o1 -o2 [options]";
  cerr << endl << endl <<
	"-a, --input1=STRING             Input filename for metagenome group A" << endl <<
	"-o, --output1=STRING            Prefix for output of metagenome A" << endl <<
  "-I, --cluster_iteration=INT           number of iteration for LSH <default 100>" << endl <<
  "-N, --min_similarity=FLOAT           minimum threshold of similarity <default 0.80>" << endl <<
  "-X, --min_similarity=FLOAT           maximum threshold of similarity <default 0.95>" << endl <<
	"-T, --threads_to_use=INT        Number of threads for running KMC etc. <default 8>" << endl <<
	"-R, --max-memory=INT            Max memory for running KMC <default 12>" << endl <<
	"    --verbose                   Print messages during run" << endl << endl <<
	;
}

void SetHyperParams(HyperParams* params) {
  (*params).cluster_iteration = 100;  //TODO: Tune this param.
  (*params).max_similarity = 0.95;
  (*params).min_similarity = 0.80;  //TODO: Tune this param.
  (*params).threads_to_use = 12;
  (*params).max_memory = 12;
}

void ParsingCommands(int argc, char*argv[], HyperParams* params) {
  int verbose_flag = 0;
  int only_flag = 0;
  string mode = "";
  const char* opt_string = "o:i:I:N:X:R:T:";
  static struct option long_options[] =
  {
    {"verbose", no_argument,  &verbose_flag, 1},
    {"output", required_argument, 0, 'o'},
  	{"input", required_argument, 0, 'i'},
    {"cluster_iteration", optional_argument, 0, 'I'},
    {"min_similarity", optional_argument, 0, 'N'},
    {"max_similarity", optional_argument, 0, 'X'},
  	{"max-memory", optional_argument, 0, 'R'},
  	{"threads_to_use", optional_argument, 0, 'T'},
  	{0,0,0,0}
  };

    int option_index = 0;
    int c;
    stringstream ss;
    while (true) {
      c = getopt_long(argc,argv,opt_string, long_options, &option_index);

      if (c == -1) {
        break;
      }
      switch (c) {
        case 'o':
          (*params).output = optarg;
          break;
        case 'i':
          (*params).input = optarg;
  	      break;
        case 'I':
          (*params).cluster_iteration = atoi(optarg);
          break;
      	case 'N':
          (*params).min_similarity = atof(optarg);
          break;
        case 'X':
          (*params).max_similarity = atof(optarg);
          break;
        case 'R':
          (*params).max_memory = atoi(optarg);
          break;
        case 'T':
          (*params).threads_to_use = atoi(optarg);
          break;
        default:
          break;
      }
    }
    if (verbose_flag) {
      (*params).verbose = true;
    }

  }


void kmerCluster(HyperParams& params){
  vector<string> samples, kmc_names;
  int tot_sample;
  int bucket_size_threshold = 1000000;
  const uint64_t batch_thresh = 100000000;
  uint64_t extracted_size_thresh = 10000000;

  auto start_time_total = chrono::high_resolution_clock::now();

  GetInput(params.input, samples, kmc_names);


	tot_sample = samples.size();
	if (params.verbose) {
		cout << endl << "# samples : " << num_sample << endl;
	}


  if (params.clustering){
    //store v_kmer without buildKHtable
    if(!params.bin){
      v_kmers.reserve(tot_sample);
      ifstream logStream("kmer_count.log");
      string line;
      getline(logStream, line);
      istringstream ss(line);
      ss >> kmap_size;
      for (int i = 0; i < tot_sample; i++) {
        ss >> kmer_coverage;
        v_kmers.push_back(kmer_coverage/kmap_size);
      }
	  }

    vector<Abundance*> unknown_abundance, *unknown_abundance_ptr;
    unknown_abundance_ptr = & unknown_abundance;
  
    Cluster(unknown_abundance_ptr, params.min_similarity, params.max_similarity, params.cluster_iteration, params.threads_to_use,tot_sample, bucket_size_threshold, params.verbose);

    
    // Save clusters.
    auto start_time_save_result = chrono::high_resolution_clock::now();
    if(params.verbose){
      cout << "Saving cluster results starts: " << endl;
    }
    IOMat::SaveResult(unknown_abundance_ptr,  params.clust_file_name+".clust", true, 5, params.verbose);
    IOMat::SaveBinary(unknown_abundance_ptr, params.clust_file_name, true, 5, params.verbose);
    auto end_time = chrono::high_resolution_clock::now();
    auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time_save_result).count();
    if(params.verbose){
      cout << "Save cluster results takes secs: " << elapsed_read << endl;
    }
    // Release allocated memory for spectra w/ charge; save pointers to spectra
    // w/o charge to variable 'spectra_of_no_charge'.
    auto start_time = chrono::high_resolution_clock::now();
    if(params.verbose){
      cout << "Releasing memory starts." << endl;
    }
    for(size_t i = 0; i< unknown_abundance.size(); i++){
        delete unknown_abundance[i];
    }

    end_time = chrono::high_resolution_clock::now();
    elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time).count();
    if(params.verbose){
      cout << "Releasing memory takes: " << elapsed_read << endl;
    }

  }
  
}

int main(int argc, char **argv) {
	if (argc < 2) {
		kLSH_PrintUsage();
	} else {
    	HyperParams params;
    	SetHyperParams(&params);
    	ParsingCommands(argc, argv, &params);
    	PrintHyperParams(params);
    	kmerCluster(params);
	}

}
