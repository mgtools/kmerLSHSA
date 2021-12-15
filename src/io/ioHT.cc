#include "ioHT.h"
//added by Mingjie on 09/05/2014
void GetInput(string input, vector<string> &samples, vector<string> &kmc_names) {
	ifstream infile (input.c_str());
	string line;
	if (infile.is_open()) {
		while (std::getline(infile, line)) {
			istringstream ss(line);
			string sample, kmc_name;
			ss >> sample >> kmc_name;

			samples.push_back(sample);
			kmc_names.push_back(kmc_name);
			//cout << sample << "\t" << kmc_name << endl;
		}
		infile.close();
	}
	else cerr << "Unable to open info file";
}


//initialize hash table with 0
void InitializeHT(ckhmap_t *kmap_ptr) {
	//int i = 0;
	for(auto it=kmap_ptr->begin(); !it.is_end(); ++it) {
		it.set_value(0);
	}
}

//added by Mingjie on 02/11/2015
void WriteHT(ckhmap_t *kmap_ptr, FILE *outfile) {
	//uint32_t kcount;
	uint16_t kcount;
	unsigned int buf_size = 256000;
	uint16_t tmp[buf_size];
	size_t i = 0;
	for(auto it=kmap_ptr->cbegin(); !it.is_end();++it) {
		kcount = it->second;
		//counter_len = Int2PChar((uint64_t)kcount, (uchar*)str);
		//fwrite((const void*) & kcount, sizeof(uint32_t), 1, outfile);
		if (i == buf_size) {
			fwrite(tmp, sizeof(uint16_t), buf_size, outfile);
			i = 0;
		}
		tmp[i] = kcount;
		//fwrite((const void*) & kcount, sizeof(uint16_t), 1, outfile);
		i ++;
	}
	if (i > 0 && i <= buf_size) {
		fwrite(tmp, sizeof(uint16_t), i, outfile);
	}

	fflush(outfile);

}

//added by Mingjie on 02/12/2015
//load kmer count info in batches
void ReadHT(std::ifstream &infile, int num_sample, uint64_t num_kmer,  uint16_t **ary_count, uint64_t batch_size, streamoff batch_offset) {
	streamoff sample_offset = 0;
	streamoff tot_offset;
	//uint32_t kcount;
	//uint16_t kcount;
	for(int i=0; i<num_sample; i++) {
		sample_offset = i*num_kmer*sizeof(uint16_t);
		tot_offset = sample_offset + batch_offset*sizeof(uint16_t);
		if(infile.is_open()) {
			infile.seekg(tot_offset, ios::beg);
			infile.read(reinterpret_cast<char*> (&ary_count[i][0]), sizeof(ary_count[0][0])*batch_size);
			//cout << ary_count[i][0] << '\t' << ary_count[i][batch_size-1] << endl;
			/*for(uint64_t j=0;j<batch_size;j++){
				//infile.read(reinterpret_cast<char*>(&kcount), sizeof(uint32_t));
				infile.read(reinterpret_cast<char*>(&kcount), sizeof(uint16_t));
				vec_count[i][j] = kcount;
			}*/
		}
		//infile.seekg(0);
		//infile.clear();
	}
	//exit(1);
}

void buildKHtable(vector<float_t>* v_kmers, size_t* kmap_size,  bool kmc, bool verbose, size_t ksize, int count_min, unsigned int num_threads, int max_memory, vector<string> samples, vector<string> kmc_names){
//void buildKHtable(vector<float_t>* v_kmers, size_t* kmap_size, pool &tp, bool kmc, bool verbose, size_t ksize, int count_min, unsigned int num_threads, int max_memory, vector<string> samples, vector<string> kmc_names){
	int tot_sample = samples.size();
	float_t kmer_coverage;
	ckhmap_t kmap, *kmap_ptr;
	kmap_ptr = &kmap;
	auto start_time_total = chrono::high_resolution_clock::now();

	//call kmc
	if (kmc) {
		if (verbose) {
			cout << endl << "...Running KMC using " << num_threads << " threads..." << endl;
      cout << "max_memory: " << max_memory << endl;
      cout << "kmer length: " << ksize << endl;
    }
    auto start_time = chrono::high_resolution_clock::now();

    for (int i=0; i<tot_sample; i++) {
      string cmd = "kmc -k" + std::to_string(ksize) + " -r -cs65535 -ci"+std::to_string(count_min) +" -t" + std::to_string(num_threads) + " -m" + std::to_string(max_memory) + " " + samples[i] + " " + kmc_names[i] + " .";
      std::system(cmd.c_str());
    }
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::minutes>(end_time - start_time).count();
    cout << "Running KMC takes " << duration << " minutes\n";
  }
  //exit(0);

  auto t1 = chrono::high_resolution_clock::now();

  if (verbose) {
    cout << endl << "...Reading kmc files..." << endl;
  }

  for (int i=0; i<tot_sample; i++) {
	  cout << kmc_names[i] << " : start reading" << endl;
	  KmcRead(kmc_names[i], kmap_ptr, verbose, num_threads, ksize);
    //KmcRead(kmc_names[i], kmap_ptr, verbose, tp, num_threads);
	  cout << "finishing reading kmc" << endl;
    if (verbose) {
      cout << kmap.size() << endl;
      cout << "load factor: " << kmap.load_factor() << endl;
    }
  }

  //record the kmer in vector based on the order of kmer in the hash table
  //for random access of kmers in BatchTest

  auto t2 = chrono::high_resolution_clock::now();
  auto duration12 = duration_cast<std::chrono::minutes>(t2 - t1).count();
  cout << "Running KmcRead takes " << duration12 << " minutes\n";

  FILE *kmer_file;
  if ((kmer_file = fopen("kmer_set.hex", "w")) == NULL) {
    cerr << "Unable to write kmer_set.hex file\n";
    exit(1);
  }

  (*kmap_size) = kmap.size();

  Kmer km;
  size_t idx_k = 0;
  for(auto it=kmap.cbegin(); !it.is_end();++it) {
    km = it->first;
    km.writeBytes(kmer_file);
    idx_k ++;
  }
  fclose(kmer_file);

  auto t3 = chrono::high_resolution_clock::now();
  auto duration23 = duration_cast<std::chrono::minutes>(t3 - t2).count();
  cout << "Creating vector and writing kmer_set.hex takes " << duration23 << " minutes\n";

  //write out kmer count info
  FILE *bin_count_file;
  FILE *log_file; //record total number of kmers in each sample
  if ((bin_count_file = fopen("kmer_count.bin", "wb")) == NULL) {
    cerr << "Unable to write kmer_count.bin file\n";
    exit(1);
  }
  if ((log_file = fopen("kmer_count.log", "w")) == NULL) {
    cerr << "Unable to write kmer_count.log file\n";
    exit(1);
  }

  if (verbose) {
    cout << "\n...Counting kmers and writing to binary..." << endl;
  }

  fprintf(log_file, "%llu", (*kmap_size));

  for(vector<string>::iterator it=kmc_names.begin(); it!=kmc_names.end(); ++it) {
    string tmp = *it;

    if (verbose) {
			cout << tmp << endl;
    }

    InitializeHT(kmap_ptr);

    kmer_coverage = KmcCount(tmp, kmap_ptr, verbose, num_threads, ksize);
    //kmer_coverage = KmcCount(tmp, kmap_ptr, verbose, tp, num_threads);
    v_kmers->push_back(kmer_coverage/(*kmap_size));
    fprintf(log_file, "\t%f", kmer_coverage);

    WriteHT(kmap_ptr, bin_count_file);
  }
  fclose(bin_count_file);
  fclose(log_file);

  auto t4 = chrono::high_resolution_clock::now();
  auto duration34 = duration_cast<std::chrono::minutes>(t4 - t3).count();
  if(verbose){
    cout << "Running KmcCount and writing to binary takes " << duration34 << " minutes\n";
  }
  //free the memory of kmap
  kmap.clear();
}
