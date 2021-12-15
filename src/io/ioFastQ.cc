#include "ioFastQ.h"

namespace Utility {

  void IOFQ::CheckRead(uset_t *g_kmer_ptr, vector<ReadEntry> &read_vec, vector<int> &record_vec, unsigned int num_threads, int tid, float kmer_vote) {

  	size_t k = Kmer::k;
  	size_t num_reads = record_vec.size();
  	float kmer_count;
  	float factor;
  	float tmp_ratio;
  	float factor_ratio;
  	float random;
  	uset_t::iterator it;
  	//uset_t::iterator it_k;

  	for (size_t i = 0; i < num_reads; i++) {
  		if ( (i % num_threads) == tid ) {
  			ReadEntry read_entry = read_vec[i];
  			if (*(read_entry.s) == '\0') {
  				cout << "\nabnormal read entry skipped\n";
  				cout << read_entry.name << endl;
  				continue;
  			}
  			if (read_entry.len < k+10) continue; //at leaset 10 kmers
  			//cout << read_entry.name << '\t' << read_entry.len << endl;
  			//cout << read_entry.s <<endl;
  			Kmer km(read_entry.s);
  			kmer_count = 0;
  			factor = 0.5;
  			for (size_t j = 0; j <= (read_entry.len - k); j++) {
  				if (j > 0) {
  					km = km.forwardBase((read_entry.s)[j+k-1]);
  				}
  				Kmer tw = km.twin();
  				Kmer rep = (km < tw) ? km : tw;

  				//char tmp[1024];
  				//rep.toString(tmp);
  				//cout << tmp << endl;
  				it = g_kmer_ptr->find(rep);
  				if (it != g_kmer_ptr->end()) {
            		kmer_count++;
          /*  it_k = large_kmer_count_ptr->find(rep);
            if (it_k != large_kmer_count_ptr->end()){
              factor = factor + 1;
              //kmer_count = kmer_count + 0;
            }
            else{
              kmer_count ++;
            }
            */
  					//char tmp[1024];
  					//km.toString(tmp);
  					//cout << tmp << endl;
  				}
  			}
  			//exit(1);
  			tmp_ratio = kmer_count/(read_entry.len - k+1);
			//cout << "tmp_ratio  : "<< tmp_ratio << endl;

  			//cout << kmer_count << endl << read_entry.len << endl << k << endl;
  			//cout << tmp_ratio << endl;
  			//exit(1);
  			if (tmp_ratio > kmer_vote) {
          		record_vec[i] = 1;
          /* randomly extract
  				factor_ratio = factor/kmer_count;
  				random = ((float) rand() / (RAND_MAX));
  				//cout << tmp_ratio << endl;
  				if(random > factor_ratio){
  					record_vec[i] = 1;
  				}*/
  			}
  		}
  	}
  }
  //extract distinctive reads for each sample
  void IOFQ::ReadExtract(uset_t *g_kmer_ptr, vector<string> &files, string output, float kmer_vote, bool verbose, unsigned int num_threads) {
  //void IOFQ::ReadExtract(uset_t *g_kmer_ptr, vector<string> &files, string output, float kmer_vote, bool verbose, pool &tp, unsigned int num_threads) {
  	//size_t num_reads;
  	FastqFile FQ(files);
  	FILE *of = fopen(output.c_str(), "w");
  	if (of == NULL) {
  		cerr << "Could not open file for writing, " << output << endl;
  		exit(1);
  	}

  	/*pool tp(num_threads);
  	if (verbose) {
  		cout << "start " << num_threads << " threads" << endl;
  	}*/

  	vector<ReadEntry> read_vec(FQ.part_size);
  	size_t part_num_reads;

  	while(FQ.read(read_vec, part_num_reads) > 0) {
  		vector<int> check_vec(part_num_reads, 0);

		#pragma omp parallel for
  		for (unsigned int tid = 0; tid < num_threads; tid++) {
  			CheckRead( g_kmer_ptr, read_vec, check_vec, num_threads, tid, kmer_vote);
			//tp.schedule(boost::bind(CheckRead, g_kmer_ptr, boost::ref(read_vec), boost::ref(check_vec), num_threads, tid, kmer_vote));
  		}

  		//tp.wait(); // wait until all kmers in the read are checked

  		//for (int i = 0; i < part_num_reads; i++) {
  		//	cout << check_vec[i] << "  " << endl;
  		//}
  		//exit(1);
  		int buf_size_factor = 100;
  		char str[8192*buf_size_factor];
  		int tot_write_len = 0;
  		int start_idx = 0;
  		int qual_len;
  		for (size_t i = 0; i < part_num_reads; i++) {
  			if (check_vec[i] == 1) {
  				ReadEntry short_read = read_vec[i];
  				//fprintf(of, "@%s\n", short_read.name);
  				//fprintf(of, "%s\n", short_read.s);
  				//fprintf(of, "+\n");
  				//fprintf(of, "%s\n", short_read.qual);
  				str[start_idx] = '@';
  				memcpy(&str[start_idx + 1], short_read.name, (short_read.name_len)*sizeof(char));
  				str[start_idx + short_read.name_len + 1] = '\n';

  				memcpy(&str[start_idx + short_read.name_len + 2], short_read.s, (short_read.len)*sizeof(char));
  				str[start_idx + short_read.name_len + short_read.len + 2] = '\n';

  				memcpy(&str[start_idx + short_read.name_len + short_read.len + 3], "+\n", 2);

  				qual_len = strlen(short_read.qual);
  				memcpy(&str[start_idx + short_read.name_len + short_read.len + 5], short_read.qual, qual_len*sizeof(char));
  				str[start_idx + short_read.name_len + short_read.len + qual_len + 5] = '\n';

  				tot_write_len += (short_read.name_len + short_read.len + qual_len + 6);

  				start_idx = tot_write_len;
  				if (tot_write_len > 8192*(buf_size_factor-1)) {
  					//cout << "here\n";
  					fwrite(&str, 1, tot_write_len, of);
  					tot_write_len = 0;
  					start_idx = 0;
				}

  				//fwrite(&str, 1, short_read.name_len + short_read.len + qual_len + 5, of);
			}	
  		}
  		if (tot_write_len > 0 && tot_write_len <= 8192*(buf_size_factor-1))
  			fwrite(&str, 1, tot_write_len, of);
  		check_vec.clear();
  		read_vec.clear();
  		if (part_num_reads < FQ.part_size) break;

  	}
  	fclose(of);
  	FQ.close();
  }
  
  
  void IOFQ::Extracting(vector<string> samples, uset_t *g_kmer_ptr, string out, int num_threads, float kmer_vote, bool verbose){
	string filename;
	auto start_time = chrono::high_resolution_clock::now();
    //unsigned int write_threads;
    //write_threads = opt.num_threads > cap_threads ? opt.num_threads : act_threads;
    //pool tp(num_threads);

    if (verbose) {
      cout << "start " << num_threads << " threads" << endl;
    }

    for(vector<string>::iterator it=samples.begin(); it!=samples.end(); ++it) {
      vector<string> tmp;
      tmp.push_back(*it);

      const char * basep = (*it).c_str();
      std::string base(basename(basep));


      filename = std::string(out) + ("_" + base);
      if (verbose) {
        cout << "writing to " << filename << endl;
      }
      ReadExtract(g_kmer_ptr, tmp, filename, kmer_vote, verbose, num_threads);
	  //ReadExtract(g_kmer_ptr, tmp, filename, kmer_vote, verbose, tp, num_threads);
      tmp.clear();
    }

    auto end_time = chrono::high_resolution_clock::now();
    auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time).count();
	if(verbose){
    	cout << "extracting reads takes (secs): " << elapsed_read << endl;
	}
	//tp.clear();
}


}  // namespace Utility
