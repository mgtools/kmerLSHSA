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

#include "../common/Common.hpp"
//#include "utils.hpp"
#include "../utils/fastq.hpp"
#include "../kmer/Kmer.hpp"
#include "../hash/HashTables.hpp"


#include "../utils/alglib-3.15.0/src/statistics.h"
#include "../utils/alglib-3.15.0/src/ap.h"

#include "../kmer/kmc_api/kmc_file.h"
#include "../kmer/kmc_reader.hpp"

#include "../utils/threadpool.hpp"
#include <boost/bind.hpp>
//#include <boost/thread.hpp>
using namespace std;
//using namespace boost;
using namespace boost::threadpool;
using namespace std::chrono;

// structs for getopt
struct SA_ProgramOptions {
  size_t k;
  string input1;
  string input2;
  string output1;
  string output2;
  int mode;
  unsigned int num_threads;
  int max_memory;
  int count_min;
  bool verbose;
  bool kmc;
  bool bin;
  float kmer_vote; //ratio of unique kmers in a read to be claimed as unique read
  float pvalue_thresh;
  vector<string> files;

  SA_ProgramOptions() : k(0), mode(1), num_threads(8), max_memory(12), count_min(2), verbose(false), kmc(false), bin(true), kmer_vote(0.5), pvalue_thresh(0.005) {}
};

void SA_PrintUsage() {
  cerr << "SA " << SA_VERSION << endl << endl;
  cerr << "Compare frequencies of k-mers [from KMC library] and saves differential reads results for each metagenome in comparison" << endl << endl;
  cerr << "Usage: SA -k -a -b -o -p [options]";
  cerr << endl << endl <<
	"-k, --kmer-size=INT             Size of k-mers, at most " << (int) (Kmer::MAX_K-1)<< endl <<
	"-a, --input1=STRING             Input filename for metagenome group A" << endl <<
	"-b, --input2=STRING             Input filename for metagenome group B" << endl <<
	"-o, --output1=STRING            Prefix for output of metagenome A" << endl <<
	"-p, --output2=STRING            Prefix for output of metagenome B" << endl <<
	"-v, --kmer-vote=FLOAT           Percentage threshold of differential k-mers in distinctive reads <default 0.5>" << endl <<
	"-t, --number-threads=INT        Number of threads for running KMC etc. <default 8>" << endl <<
	"-m, --max-memory=INT            Max memory for running KMC <default 12>" << endl <<
	"-c, --count-min=INT            Min threshold of k-mer count for running KMC <default 2>" << endl <<
	"-q, --pvalue-thresh=FLOAT       For U-test <default 0.005>" << endl <<
	"-M, --mode=INT                  1 or 2 <default 1>" << endl <<
	"    --kmc                       Call KMC (if you do not have kmc output; make sure kmc is available in the current environment)" << endl <<
	"    --no-bin                    Do not generate kmer_count.bin (already existed)" << endl <<
    "    --verbose                   Print messages during run" << endl << endl
    ;
}


void SA_ParseOptions(int argc, char **argv, SA_ProgramOptions &opt) {
  int verbose_flag = 0;
  int kmc_flag = 0;
  int bin_flag = 0;
  const char* opt_string = "k:o:p:a:b:v:q:m:c:t:M:";
  static struct option long_options[] =
  {
    {"verbose", no_argument,  &verbose_flag, 1},
	{"kmc", no_argument, &kmc_flag, 1},
	{"no-bin", no_argument, &bin_flag, 1},
	{"kmer-size", required_argument, 0, 'k'},
    {"output1", required_argument, 0, 'o'},
	{"output2", required_argument, 0, 'p'},
	{"input1", required_argument, 0, 'a'},
	{"input2", required_argument, 0, 'b'},
	{"kmer-vote", required_argument, 0, 'v'},
	{"pvalue-thresh", required_argument, 0, 'q'},
	{"max-memory", required_argument, 0, 'm'},
	{"count-min", required_argument, 0, 'c'},
	{"num-threads", required_argument, 0, 't'},
    {"mode", required_argument, 0, 'M'},
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
	case 'k':
	  opt.k = atoi(optarg);
	  break;
	case 'o':
      opt.output1 = optarg;
      break;
	case 'p':
	  opt.output2 = optarg;
	  break;
	case 'a':
	  opt.input1 = optarg;
	  break;
	case 'b':
	  opt.input2 = optarg;
	  break;
	case 'v':
	  opt.kmer_vote = atof(optarg);
	  break;
	case 'q':
	  opt.pvalue_thresh = atof(optarg);
	  break;
	case 'M':
	  opt.mode = atoi(optarg);
	  break;
    	case 't':
	  opt.num_threads = atoi(optarg);
	  break;
	case 'm':
	  opt.max_memory = atoi(optarg);
	  break;
	case 'c':
          opt.count_min = atoi(optarg);
          break;
	default: break;
    }
  }

  if (verbose_flag) {
	  opt.verbose = true;
  }

  if (kmc_flag) {
	  opt.kmc = true;
  }

  if (bin_flag) {
	  opt.bin = false;
  }
}


bool SA_CheckOptions(SA_ProgramOptions &opt) {
  bool ret = true;

  if (opt.k <= 0 || opt.k >= MAX_KMER_SIZE) {
	  cerr << "Error, invalid value for kmer-size: " << opt.k << endl;
	  cerr << "Values must be between 1 and " << (MAX_KMER_SIZE-1) << endl;
	  ret = false;
  }


  /*
  if (opt.files1.size() == 0 || opt.files2.size() == 0) {
    cerr << "Need to specify files for input" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    vector<string>::const_iterator it;
    int intStat;
    for(it = opt.files1.begin(); it != opt.files1.end(); ++it) {
      intStat = stat(it->c_str(), &stFileInfo);
      if (intStat != 0) {
	cerr << "Error: file not found, " << *it << endl;
	ret = false;
      }
    }
	for(it = opt.files2.begin(); it != opt.files2.end(); ++it) {
      intStat = stat(it->c_str(), &stFileInfo);
      if (intStat != 0) {
	cerr << "Error: file not found, " << *it << endl;
	ret = false;
      }
    }

  }
  */



  //TODO: check if we have permission to write to outputfile

  return ret;

}

void SA_PrintSummary(const SA_ProgramOptions &opt) {
  cout << "Mode: " << opt.mode << endl;
  cout << "kmer size: " << Kmer::k << endl;
  cout << "Using kmer_vote: " << opt.kmer_vote << endl;
  cout << "Using pvalue_thresh: " << opt.pvalue_thresh << endl;
  cout << "Writing kmer_count.bin? " << std::boolalpha << opt.bin << endl;
}

//added by Mingjie on 02/12/2015
//test kmer count in batches
void BatchTest(uint16_t **ary_count, int tot_sample, uint64_t batch_size, streamoff batch_offset, int num_sample1, int num_sample2, vector<uint64_t> &v_kmers, float pvalue_thresh, const vector<Kmer> &kvec, uset_t *g_kmer1_ptr, uset_t *g_kmer2_ptr, uset_t *large_kmer_count_ptr) {
	double freq_array1[num_sample1], freq_array2[num_sample2];
	double bothtails, lefttail, righttail, freq;
	int numZero;
	uint16_t cnt;
	uint64_t tot_cnt;
	uint64_t num_kmers;
	alglib::real_1d_array array_s1, array_s2;

	bothtails = 0;
	lefttail  = 0;
	righttail = 0;

	Kmer km;
	for(uint64_t i=0; i<batch_size; i++) {
		km = kvec[batch_offset+i];
		tot_cnt = 0;
		numZero = 0;
		for(int j=0; j<tot_sample; j++){
			cnt = ary_count[j][i];
			//cout << cnt << endl;
			tot_cnt += cnt;
			if (cnt == 0){
				numZero++;
			}
			num_kmers = v_kmers[j];
			freq = double(cnt)*1000000/num_kmers;

			if(j<num_sample1) {
				freq_array1[j] = freq;
			} else {
				freq_array2[j-num_sample1] = freq;
			}
		}
		//char tmp[1024];
		//km.toString(tmp);
		//cout << tmp << endl;
		//exit(1);
		//if (tot_cnt > tot_sample) {
		array_s1.setcontent(num_sample1, freq_array1);
		array_s2.setcontent(num_sample2, freq_array2);
		alglib::mannwhitneyutest(array_s1, num_sample1, array_s2, num_sample2, bothtails, lefttail, righttail);
		if(numZero < tot_sample * 0.7 ){
			if(lefttail <= pvalue_thresh) {
				g_kmer2_ptr->insert(km);
			} else if (righttail <= pvalue_thresh) {
				g_kmer1_ptr->insert(km);
			}
			if (tot_cnt > tot_sample*2  ){
				large_kmer_count_ptr->insert(km);
			}
		}
	}
}


void SA(int argc, char **argv) {

	SA_ProgramOptions opt;
	SA_ParseOptions(argc,argv,opt);
	ckhmap_t kmap, *kmap_ptr;
	kmap_ptr = &kmap;
  uset_t large_kmer_count, *large_kmer_count_ptr;
  large_kmer_count_ptr = &large_kmer_count;
	size_t kmer_coverage;
	size_t kmap_size;
	std::string filename;
	vector<string> samples1, samples2, samples;
	vector<string> kmc_names1, kmc_names2, kmc_names;
	vector<size_t> v_kmers; //kmer_cov for each sample
	int num_sample1, num_sample2, tot_sample;
	uset_t g_kmer1, g_kmer2, *g_kmer1_ptr, *g_kmer2_ptr;  //differential kmers for two groups in comparison
	g_kmer1_ptr = &g_kmer1;
	g_kmer2_ptr = &g_kmer2;
	vector<Kmer> kvec;

	alglib::real_1d_array array_s1, array_s2;

	int cap_threads = 6;
	steady_clock::time_point t1, t2, t3, t4, t5, t6;

	if (argc < 2) {
		SA_PrintUsage();
		exit(1);
	}

	/*
	if (!SA_CheckOptions(opt)) {
	CountBF_PrintUsage();
	exit(1);
	} */


	// set static global k-value
	Kmer::set_k(opt.k);

	if (opt.verbose) {
		SA_PrintSummary(opt);
	}

	int act_threads;
	act_threads = opt.num_threads > cap_threads ? cap_threads:opt.num_threads;

	GetInput(opt.input1, samples1, kmc_names1);
	GetInput(opt.input2, samples2, kmc_names2);

	samples = samples1;
	samples.insert(samples.end(), samples2.begin(), samples2.end());
	kmc_names = kmc_names1;
	kmc_names.insert(kmc_names.end(), kmc_names2.begin(), kmc_names2.end());

	num_sample1 = samples1.size();
	num_sample2 = samples2.size();
	tot_sample = samples.size();
	if (opt.verbose) {
		cout << endl << "# samples in group 1: " << num_sample1 << endl << "# samples in group 2: " << num_sample2 << endl;
	}
  if (opt.bin){
		pool tp(opt.num_threads);
		buildKHtable();

  }
	else {
		//read in total number of kmers in each sample
		ifstream logStream("kmer_count.log");
		string line;
		getline(logStream, line);
		istringstream ss(line);
		ss >> kmap_size;
		for (int i = 0; i < tot_sample; i++) {
			ss >> kmer_coverage;
			v_kmers.push_back(kmer_coverage);
		}
	}

	//read kmers in kvec to restore the original kmer order in kmap
	ifstream kmer_file("kmer_set.hex");
	kvec.resize(kmap_size);
	uint8_t bytes[(Kmer::MAX_K)/4];
	size_t idx_k = 0;
	for (size_t i = 0; i < kmap_size; i++) {
		kmer_file.read(reinterpret_cast<char*> (&bytes[0]), sizeof(bytes[0])*(Kmer::MAX_K)/4);
		Kmer km(bytes);
		kvec[idx_k] = km;
		idx_k ++;
	}
	t4 = steady_clock::now();

	ifstream inStream("kmer_count.bin", ios::binary);

	const uint64_t batch_thresh = 10000000; //if changed, remember to change 2D array parameter vec_count in readHT and batchTest
	uint64_t batch_size;

	//initialize 2D array
	uint16_t** ary_count = new uint16_t*[tot_sample];
	for (int i=0; i < tot_sample; i++) {
		ary_count[i] = new uint16_t[batch_thresh];
	}

	if (opt.verbose) {
	  cout << "\n...Loading and testing in batches..." << endl;
	}

	//load and test kmer count info in batches
	uint64_t kcnt_rem = kmap_size;
	streamoff batch_offset = 0;
	while(kcnt_rem >= batch_thresh) {
	  batch_size = batch_thresh;

	  ReadHT(inStream, tot_sample, kmap_size, ary_count, batch_size, batch_offset);
	  BatchTest(ary_count, tot_sample, batch_size, batch_offset, num_sample1, num_sample2, v_kmers, opt.pvalue_thresh, kvec, g_kmer1_ptr, g_kmer2_ptr, large_kmer_count_ptr);

	  //cout << g_kmer1.size() << "\t" << g_kmer2.size() <<endl;

	  batch_offset += batch_size;
	  kcnt_rem -= batch_size;
	  if (opt.verbose) {
		  cout << "# loaded kmers: " << batch_offset << endl;
	  }


	}
	if(kcnt_rem > 0 && kcnt_rem < batch_thresh) {
	  batch_size = kcnt_rem;
	  ReadHT(inStream, tot_sample, kmap_size, ary_count, batch_size, batch_offset);
	  BatchTest(ary_count, tot_sample, batch_size, batch_offset, num_sample1, num_sample2, v_kmers, opt.pvalue_thresh, kvec, g_kmer1_ptr, g_kmer2_ptr, large_kmer_count_ptr);
	}

	for (int i = 0; i < tot_sample; ++i) {
		delete [] ary_count[i];
	}
	delete [] ary_count;
	inStream.close();

	t5 = steady_clock::now();
	auto duration45 = duration_cast<std::chrono::minutes>(t5 - t4).count();
	cout << "BatchTesting takes " << duration45 << " minutes\n";

	if (opt.verbose) {
	  cout << "\n# differential kmers in group 1: " << g_kmer1.size() << endl;
	  cout << "# differential kmers in group 2: " << g_kmer2.size() << endl;
	}

	//deallocate the memory of kvec
	vector<Kmer>().swap(kvec);

	if (opt.verbose) {
	  cout << "\n...Writing distinctive reads..."  << endl;
	}


	//unsigned int write_threads;
	//write_threads = opt.num_threads > cap_threads ? opt.num_threads : act_threads;
	pool tp(opt.num_threads);

	if (opt.verbose) {
		cout << "start " << tp.size() << " threads" << endl;
	}

	for(vector<string>::iterator it=samples1.begin(); it!=samples1.end(); ++it) {
	  vector<string> tmp;
	  tmp.push_back(*it);

	  const char * basep = (*it).c_str();
	  std::string base(basename(basep));


		filename = std::string(opt.output1) + ("_" + base);
		if (opt.verbose) {
			cout << "writing to " << filename << endl;
		}
		ReadExtract(g_kmer1_ptr, tmp, filename, opt.kmer_vote, opt.verbose, tp, opt.num_threads, large_kmer_count_ptr);
	  tmp.clear();
	}

	for(vector<string>::iterator it=samples2.begin(); it!=samples2.end(); ++it) {
	  vector<string> tmp;
	  tmp.push_back(*it);

	  const char * basep = (*it).c_str();
	  std::string base(basename(basep));
		filename = std::string(opt.output2) + ("_" + base);
		if (opt.verbose) {
			cout << "writing to " << filename << endl;
		}
		ReadExtract(g_kmer2_ptr, tmp, filename, opt.kmer_vote, opt.verbose, tp, opt.num_threads, large_kmer_count_ptr);
	  tmp.clear();
	}

	t6 = steady_clock::now();
	auto duration56 = duration_cast<std::chrono::minutes>(t6 - t5).count();
	cout << "ReadExtract takes " << duration56 << " minutes\n";

}
