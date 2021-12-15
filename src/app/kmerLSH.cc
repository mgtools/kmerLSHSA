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
#include "../io/ioHT.h"
#include "../io/ioFastQ.h"

#include "../utils/fastq.h"

#include "../hash/HashTables.h"
#include "../hash/hash.h"
#include "../hash/lshash.h"

#include "../kmer/kmc_api/kmc_file.h"
#include "../kmer/kmc_reader.h"
#include "../kmer/Kmer.h"

//#include "../utils/threadpool.hpp"
//#include <boost/bind.hpp>

#include "../utils/alglib-3.17.0/src/statistics.h"
#include "../utils/alglib-3.17.0/src/ap.h"
//#include <boost/thread.hpp>

using namespace std;
//using namespace boost;
//using namespace boost::threadpool;
using namespace std::chrono;
using namespace Core;
using namespace Utility;

struct HyperParams
{
  size_t k;
  float min_similarity, max_similarity;
  float scale;
  int cluster_iteration;
  int max_memory;
  int count_min;
  int size_thresh;
  uint64_t batch_size, extracted_size;
  float pval_thresh;
  float kmer_vote;
  bool verbose;
  bool kmc, bin, clustering, extracting;
  unsigned int threads_to_use;
  string tmp_dir;
  string input1, input2;
  string output1, output2;
  string clust_file_name;
};

void PrintHyperParams(const HyperParams& params) {
  cout << "************ kmers Cluster for Subtractiva Assembly Params Setting ****************" << endl;
  cout << "cluster iteration: " << params.cluster_iteration << endl;
  cout << "group1 result prefix: " << params.output1 << endl;
  cout << "group2 result prefix: " << params.output2 << endl;
  cout << "cluster result file: " << params.clust_file_name << endl;
  cout << "group1 input file: " << params.input1 << endl;
  cout << "group2 input file: " << params.input2 << endl;
  cout << "min similarity: " << params.min_similarity << endl;
  cout << "max similarity: " << params.max_similarity << endl;
  cout << "threds to use: " << params.threads_to_use << endl;
  cout << "cluster size threshold: " << params.size_thresh << endl;
  cout << "percentage threshold of differential k-mers in distinctive reads" << params.kmer_vote << endl;
  cout << "********************************************************" << endl;
}

void kLSH_PrintUsage() {
  cerr << "kmerLSHSA "<< KLSH_VERSION << endl << endl;
  cerr << "Clustering of k-mers for Subtractive Assembly " << endl << endl;
  cerr << "Usage: ./kmerLSHSA -a input1 -b input2 -o output1_prefix -p output2_prefix [options]";
  cerr << endl << endl <<
	"-a, --input1=STRING             Input filename for metagenome group A" << endl <<
  "-b, --input2=STRING             Input filename for metagenome group B" << endl <<
	"-o, --output1=STRING            Prefix for output of metagenome A" << endl <<
  "-p, --output2=STRING            Prefix for output of metagenome B" << endl <<
  "-I, --cluster_iteration=INT           number of iteration for LSH <default 100>" << endl <<
  "-N, --min_similarity=FLOAT           minimum threshold of similarity <default 0.85>" << endl <<
  "-X, --min_similarity=FLOAT           maximum threshold of similarity <default 0.95>" << endl <<
  "-K, --kmer_size=INT             Size of k-mers, at most " << (int) (Kmer::MAX_K-1)<< endl <<
	"-T, --threads_to_use=INT        Number of threads for running KMC etc. <default 8>" << endl <<
	"-R, --max-memory=INT            Max memory for running KMC <default 12>" << endl <<
	"-C, --count-min=INT            Min threshold of k-mer count for running KMC <default 2>" << endl <<
  "-S, --size_thresh=INT       Threshold of the size of clustering for T-Test <default 500000>" << endl <<
  "-V, --kmer_vote=FLOAT           Percentage threshold of differential k-mers in distinctive reads <default 0.5>" << endl <<
  "-F, --clust_file_name=STRING           intermediate clustering result file name <default 'clustering_result.txt'>" << endl <<
	"-M, --mode=STRING                Optional K : run kmc, B : make bin file, C : clustering, E : extract differential reads " << endl <<
	"    --verbose                   Print messages during run" << endl << endl <<
	"    --only                   Run only the setting mode " << endl << endl
	;
}

void SetHyperParams(HyperParams* params) {
  (*params).cluster_iteration = 100;  //TODO: Tune this param.
  (*params).max_similarity = 0.95;
  (*params).min_similarity = 0.85;  //TODO: Tune this param.
  (*params).scale = 1;
  (*params).threads_to_use = 12;
  (*params).kmc = true;
  (*params).bin = true;
  (*params).extracting = true;
  (*params).clustering = true;
  (*params).max_memory = 12;
  (*params).count_min = 2;
  (*params).k = 23;
  (*params).size_thresh = 500000;
  (*params).kmer_vote = 0.5;
  (*params).tmp_dir = "tmp/";
  (*params).clust_file_name = "clustering_result.txt";
}

void ParsingCommands(int argc, char*argv[], HyperParams* params) {
  int verbose_flag = 0;
  int only_flag = 0;
  string mode = "";
  const char* opt_string = "o:p:a:b:H:I:N:X:C:T:K:S:P:V:F:M:";
  static struct option long_options[] =
  {
    {"verbose", no_argument,  &verbose_flag, 1},
    {"only", no_argument,  &only_flag, 1},
    {"output1", required_argument, 0, 'o'},
    {"output2", required_argument, 0, 'p'},
  	{"input1", required_argument, 0, 'a'},
    {"input2", required_argument, 0, 'b'},
    {"cluster_iteration", optional_argument, 0, 'I'},
    {"min_similarity", optional_argument, 0, 'N'},
    {"max_similarity", optional_argument, 0, 'X'},
  	{"max-memory", optional_argument, 0, 'R'},
  	{"count-min", optional_argument, 0, 'C'},
  	{"threads_to_use", optional_argument, 0, 'T'},
    {"kmer_size", optional_argument, 0, 'K'},
    {"size_thresh", optional_argument, 0, 'S'},
    {"kmer_vote", optional_argument, 0, 'V'},
    {"tmp_dir", optional_argument, 0, 'D'},
    {"clust_file_name", optional_argument, 0, 'F'},
    {"mode", optional_argument, 0, 'M'},
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
          (*params).output1 = optarg;
          break;
        case 'p':
          (*params).output2 = optarg;
          break;
        case 'a':
          (*params).input1 = optarg;
  	      break;
        case 'b':
          (*params).input2 = optarg;
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
        case 'C':
          (*params).count_min = atoi(optarg);
          break;
        case 'T':
          (*params).threads_to_use = atoi(optarg);
          break;
        case 'K':
          (*params).k = atoi(optarg);
          break;
        case 'S':
          (*params).size_thresh = atoi(optarg);
          break;
        case 'V':
          (*params).kmer_vote = atof(optarg);
          break;
        case 'D':
          (*params).tmp_dir = optarg;
          break;
        case 'F':
          (*params).clust_file_name = optarg;
          break;
        case 'M':
          mode = optarg;
          break;
        default:
          break;
      }
    }
    if (verbose_flag) {
      (*params).verbose = true;
    }

    if (only_flag){
      if (mode == "K"){
        (*params).bin = false;
        (*params).clustering = false;
        (*params).extracting = false;
      }
      else if (mode == "B"){
        (*params).kmc = false;
        (*params).clustering = false;
        (*params).extracting = false;
      }
      else if (mode == "C"){
        (*params).bin = false;
        (*params).kmc = false;
        (*params).extracting = false;
      }
      else if (mode == "E"){
        (*params).bin = false;
        (*params).kmc = false;
        (*params).clustering = false;
      }
    } else {
      if (mode == "B"){
        (*params).kmc = false;
      }
      else if (mode == "C"){
        (*params).bin = false;
        (*params).kmc = false;
      }
      else if (mode == "E"){
        (*params).bin = false;
        (*params).kmc = false;
        (*params).clustering = false;
      }
    }
  }

void init_clustering(vector<Abundance*>* unknown_abundance_ptr, int dim, float min_similarity, float max_similarity, int cluster_iteration, unsigned int threads_to_use, int tot_sample, size_t kmap_size, vector<float_t> v_kmers, int bucket_size_threshold, string tmp_dir, const uint64_t batch_thresh, bool verbose){
  auto start_time = chrono::high_resolution_clock::now();
  ifstream inStream("kmer_count.bin", ios::binary);
  vector<Abundance*> out_abundance, local_abundance, *out_abundance_ptr, *local_abundance_ptr;
  local_abundance_ptr = & local_abundance;
  out_abundance_ptr = & out_abundance;
  float similarity = min_similarity;
  //const uint64_t batch_thresh = 100000000; //if changed, remember to change 2D array parameter vec_count in readHT and batchTest
  uint64_t batch_size;
  string read_tmp_file = "", write_tmp_file = "", read_tmp_file_clust = "", write_tmp_file_clust = "";
  deque<int> tmp_files;
  size_t total_size = 0;

  //initialize 2D array
  uint16_t** ary_count = new uint16_t*[tot_sample];
  for (int i=0; i < tot_sample; i++) {
  	ary_count[i] = new uint16_t[batch_thresh];
  }

  if (verbose) {
    cout << "\n...Loading and testing in batches..." << endl;
	  cout << v_kmers.size() << endl;
  }

  //load and test kmer count info in batches
  uint64_t kcnt_rem = kmap_size;
  streamoff batch_offset = 0;
  int iter = kmap_size / batch_thresh;
  int tmp = 0;
  if(verbose){
    cout << "iteration : " << iter << " kmap_size : " << kmap_size << endl;
  }
  //First, read the binary file of k-mers hash table. 
  //Clustering the block of hash table
  //Saving to the other binary file.
  for(int i = 0; i< iter+1 ; i++){
    if(i == iter){
      //batch_offset = batch_thresh*i;
      batch_size = kcnt_rem - batch_offset;
    }
    else{
      batch_size = batch_thresh;
      //batch_offset = batch_thresh * i;
    }
    //cout << "i: " << i << " batch_size : " << batch_size << " batch_offset : " << batch_offset << endl;
    ReadHT(inStream, tot_sample, kmap_size, ary_count, batch_size, batch_offset);
    IOMat::convertHTMat(ary_count, v_kmers,  tot_sample, verbose, batch_size, batch_offset, local_abundance_ptr);
    Cluster(local_abundance_ptr, similarity, max_similarity,  cluster_iteration,  threads_to_use, dim, bucket_size_threshold, verbose);
    total_size += local_abundance.size();
    //////write result on tmp bin
    if(i == 0){
      write_tmp_file = tmp_dir+to_string(tmp)+".bin";
      write_tmp_file_clust = write_tmp_file+".clust";
      IOMat::SaveResult(local_abundance_ptr,  write_tmp_file_clust, true, 0, verbose);
      IOMat::SaveBinary(local_abundance_ptr, write_tmp_file, true, 0, verbose);
      tmp++;
    }
    else{
      IOMat::SaveResult(local_abundance_ptr,  write_tmp_file_clust, false, 0, verbose);
      IOMat::SaveBinary(local_abundance_ptr, write_tmp_file, false,  0, verbose);
    }
	  for(size_t j = 0; j< local_abundance.size(); j++){
      delete local_abundance[j];
    }
    vector<Abundance*>().swap(local_abundance);
	  batch_offset += batch_size;
    if (verbose) {
  	  cout << "# loaded kmers: " << batch_offset << endl;
  	}
  }
  
  for (int i = 0; i < tot_sample; ++i) {
  	delete [] ary_count[i];
  }
  delete [] ary_count;
  inStream.close();

  /////read temporary directory files/////////
  while (total_size > batch_thresh*5){
    similarity -=0.001;
    batch_offset = 0;
    read_tmp_file = write_tmp_file;
    read_tmp_file_clust = write_tmp_file_clust;
    write_tmp_file = tmp_dir+to_string(tmp)+".bin";
    write_tmp_file_clust = write_tmp_file+".clust";
    tmp++;
    iter = total_size / batch_thresh;
    kcnt_rem = total_size;
    total_size = 0;
    for(int i = 0; i< iter+1; i++){
      if (i == iter ){
        //batch_offset = batch_thresh*i;
        batch_size = kcnt_rem - batch_offset;
      }
      else{
        //batch_offset = batch_thresh*i;
        batch_size = batch_thresh;
      }
      cout << "i: " << i << " batch_size : " << batch_size << " batch_offset : " << batch_offset << endl;
      IOMat::ReadCluster(local_abundance_ptr, dim, read_tmp_file, batch_offset, batch_size, verbose);
      //IOMat::ReadCluster(local_abundance_ptr, read_tmp_file+".clust", batch_offset, batch_size, verbose);
      Cluster(local_abundance_ptr, similarity, max_similarity, cluster_iteration+4,  threads_to_use, dim, bucket_size_threshold, verbose);
      total_size += local_abundance.size();

      if(i == 0){
        IOMat::SaveResult(local_abundance_ptr,  write_tmp_file_clust, true, 0, verbose);
        IOMat::SaveBinary(local_abundance_ptr, write_tmp_file, true, 0, verbose);
      }
      else{
        IOMat::SaveResult(local_abundance_ptr,  write_tmp_file_clust, false, 0, verbose);
        IOMat::SaveBinary(local_abundance_ptr, write_tmp_file, false, 0, verbose);
      }
    
      for(size_t j = 0; j< local_abundance.size(); j++){
        delete local_abundance[j];
      }
      vector<Abundance*>().swap(local_abundance);
	    batch_offset += batch_size;
	    if (verbose) {
      	cout << "# loaded kmers: " << batch_offset << endl;
      }
    }
    ///////////delete read temporary file ///////////////
    if(remove(read_tmp_file.c_str()) != 0){
      perror("The temporary file deletion failed");
    }
    else{
      cout << read_tmp_file << "file are removed" << endl;
    }
    if(remove(read_tmp_file_clust.c_str()) != 0){
      perror("The temporary file deletion failed");
    }
    else{
      cout << read_tmp_file_clust << "file are removed" << endl;
    }
  }
  
  read_tmp_file = write_tmp_file;
  read_tmp_file_clust = write_tmp_file_clust;
  IOMat::ReadClusterAll(local_abundance_ptr, dim, read_tmp_file, verbose);
  //IOMat::ReadClusterAll(local_abundance_ptr, read_tmp_file+".clust");

  unknown_abundance_ptr->swap(local_abundance);
  vector<Abundance*>().swap(local_abundance);

  
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time).count();
  
  if(verbose){
    cout << "Finish conversion matrix" << endl;
	  cout << "reading Matrix takes secs:\t" << elapsed_read << endl;
  }
  
}

void kmerCluster(HyperParams& params){
  vector<string> samples1, samples2, kmc_names1, kmc_names2 ;
  vector<string> samples, kmc_names;
  vector<float_t> v_kmers;
  int tot_sample, num_sample1, num_sample2;
  Kmer::set_k(params.k);
  size_t kmap_size;
  float_t kmer_coverage;
  int bucket_size_threshold = BUCKET_SIZE_THRESH;
  const uint64_t batch_thresh = BATCH_THRESH;
  uint64_t extracted_size_thresh = EXTRACTED_SIZE_THRESH;

  auto start_time_total = chrono::high_resolution_clock::now();

  GetInput(params.input1, samples1, kmc_names1);
  GetInput(params.input2, samples2, kmc_names2);

  samples = samples1;
	samples.insert(samples.end(), samples2.begin(), samples2.end());
	kmc_names = kmc_names1;
	kmc_names.insert(kmc_names.end(), kmc_names2.begin(), kmc_names2.end());

	num_sample1 = samples1.size();
	num_sample2 = samples2.size();
	tot_sample = samples.size();
	if (params.verbose) {
		cout << endl << "# samples in group 1: " << num_sample1 << endl << "# samples in group 2: " << num_sample2 << endl;
	}

  if (params.kmc || params.bin){
     buildKHtable( &v_kmers, &kmap_size,   params.kmc, params.verbose, params.k, params.count_min, params.threads_to_use, params.max_memory, samples, kmc_names);
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
  
    // Creating a directory
    if (mkdir((params.tmp_dir).c_str(), 0777) == -1)
        cerr << "Error :  " << strerror(errno) << endl;
  
    else
        cout << "Temporary Directory created";

    vector<Abundance*> unknown_abundance, *unknown_abundance_ptr;
    unknown_abundance_ptr = & unknown_abundance;

    init_clustering(unknown_abundance_ptr, tot_sample, params.max_similarity, 0.97, 1, params.threads_to_use, tot_sample, kmap_size, v_kmers, bucket_size_threshold, params.tmp_dir, batch_thresh, params.verbose);

  
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
  if (params.extracting){
    if(params.verbose){
      cout << "Start to extract the differential reads from raw data" << endl;
    }
    auto start_time_extracting = chrono::high_resolution_clock::now();
    vector<Abundance*> clusteredab, *clusteredab_ptr;
    clusteredab_ptr = &clusteredab;

    uset_t g_kmer1, g_kmer2, *g_kmer1_ptr, *g_kmer2_ptr;  //differential kmers for two groups in comparison
	  g_kmer1_ptr = &g_kmer1;
	  g_kmer2_ptr = &g_kmer2;

    unordered_set<uint64_t> g_kmer_id1, g_kmer_id2, *g_kmer_id1_ptr, *g_kmer_id2_ptr;
    g_kmer_id1_ptr = &g_kmer_id1;
    g_kmer_id2_ptr = &g_kmer_id2;

    string head = "";
    int dim = 0;
    //read clustering result and do statistical testing(WRS)
    IOMat::ReadClusterAll(clusteredab_ptr, tot_sample, params.clust_file_name, params.verbose);
    AB::dynamicST(g_kmer_id1_ptr, g_kmer_id2_ptr, clusteredab_ptr, num_sample1, num_sample2, params.pval_thresh, params.size_thresh, extracted_size_thresh, params.verbose);
    for(size_t i; i< clusteredab.size(); i++){
      delete clusteredab[i];
    }
    
    if(params.verbose){
      cout << "# of differential kmers in group A : " << g_kmer_id1.size() << endl;
      cout << "# of differential kmers in group B : " << g_kmer_id2.size() << endl;
    }
    //read kmers in kvec to restore the original kmer order in kmap
    
    ifstream logStream("kmer_count.log");
    string line;
    getline(logStream, line);
    istringstream ss(line);
    ss >> kmap_size;
    logStream.close();

    ifstream kmer_file("kmer_set.hex");
	  //kvec.resize(kmap_size);
	  uint8_t bytes[(Kmer::MAX_K)/4];
	  size_t idx_k = 0;

	  for (size_t i = 0; i < kmap_size; i++) {
	    kmer_file.read(reinterpret_cast<char*> (&bytes[0]), sizeof(bytes[0])*(Kmer::MAX_K)/4);
	    Kmer km(bytes);
      
      if(g_kmer_id1.find(i) != g_kmer_id1.end() ){
        g_kmer1.insert(km);
      }
      else if(g_kmer_id2.find(i) != g_kmer_id2.end()){
        g_kmer2.insert(km);
      }
	    idx_k ++;
    }

    //clear the id sets and close file
    unordered_set<uint64_t>().swap(g_kmer_id1);
    unordered_set<uint64_t>().swap(g_kmer_id2);
    kmer_file.close();

    IOFQ::Extracting(samples1, g_kmer1_ptr, params.output1, params.threads_to_use, params.kmer_vote, params.verbose);
    IOFQ::Extracting(samples2, g_kmer2_ptr, params.output2, params.threads_to_use, params.kmer_vote, params.verbose);

    // Report time for all the procedures done above.
    auto end_time_extracting = chrono::high_resolution_clock::now();
    auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time_extracting - start_time_extracting).count();
    if(params.verbose){
      cout << "extracting reads takes (secs): " << elapsed_read << endl;
    }
  }

  // Report time for all the procedures done above.
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time_total).count();
  cout << "kmerLSHSA algorithm in total takes (secs): " << elapsed_read << endl;
  
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
