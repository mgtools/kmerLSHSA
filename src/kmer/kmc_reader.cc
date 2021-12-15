
#include "kmc_reader.h"
#include <cstdlib>

void InsertKmer(ckhmap_t *kmap_ptr, vector<Kmer> &kmer_vec, int num_threads, unsigned int tid) {
	size_t num_kmer = kmer_vec.size();
	Kmer km;
	//ckhmap_t::iterator it;
	
	for (size_t i = 0; i < num_kmer; i++) {
		if ( (i % num_threads == tid) ) {
			km = kmer_vec[i];

			Kmer tw = km.twin();
			Kmer rep = (km < tw) ? km : tw;
			
			kmap_ptr->insert(rep, 0);
			//kmap_ptr->insert(KmerIntPair(rep,0));

		}
	}
}
//go through a kmc library, insert kmers to a hash table 
//  and compute total coverages of kmers in the single kmc library
//void KmcRead(string kmc_file_name, ckhmap_t *kmap_ptr, bool verbose, unsigned int num_threads) {
void KmcRead(string kmc_file_name, ckhmap_t *kmap_ptr, bool verbose, unsigned int threads_to_use, size_t ksize) {
//void KmcRead(string kmc_file_name, ckhmap_t *kmap_ptr, bool verbose, pool &tp, unsigned int num_threads) {
	CKMCFile kmer_data_base;
	//int32 i;

	if (!kmer_data_base.OpenForListing(kmc_file_name)) 
		exit(EXIT_FAILURE);
	else {
		CKMCFileInfo info;

		kmer_data_base.Info(info);
		
		char str[ksize];
		uint32_t cnt;

		vector<Kmer> kmer_vec(info.total_kmers);
		
		if (verbose) {
			cout << "total kmers : "<< info.total_kmers << endl;
		}

		//TODO: add min_count_to_set and max_count_to_set
		CKmerAPI kmer_object(info.kmer_length);

		uint64 idx_k = 0;
		while (kmer_data_base.ReadNextKmer(kmer_object, cnt)) {
			kmer_object.to_string(str);
			//cout << str << endl;
			Kmer km(str);
			
			kmer_vec[idx_k] = km;
			idx_k ++;
		}
		kmer_data_base.Close();
		cout << "finish to read " << endl;
		/*pool tp(num_threads);
		if (verbose) {
			cout << "start " << num_threads << " threads" << endl;
		}*/
		
		#pragma omp parallel for num_threads(threads_to_use)
		for (unsigned int tid = 0; tid < threads_to_use; tid++) {
			InsertKmer(kmap_ptr, kmer_vec, threads_to_use, tid);
			//tp.schedule(boost::bind(InsertKmer, kmap_ptr, boost::ref(kmer_vec), num_threads, tid));
		}
		//tp.wait();
		//tp.clear();


			//Kmer tw = km.twin();
			//Kmer rep = (km < tw) ? km : tw;

			//pair<hmap_t::iterator,bool> ref = kmap_ptr->insert(KmerIntPair(rep,0));
			//kmap_ptr->insert(rep, 0);
			//kmap_ptr->insert(KmerIntPair(rep,0));
			
		//}
		
		kmer_vec.clear();
		//tp.clear();

	}
}

void SearchKmer(ckhmap_t *kmap_ptr, vector<Kmer> &kmer_vec, vector<uint32> &cnt_vec, int num_threads, unsigned int tid) {
	size_t num_kmer = kmer_vec.size();
	Kmer km;
	uint32 cnt, tmp;
	
	for (size_t i = 0; i < num_kmer; i++) {
		if ( (i % num_threads == tid) ) {
			km = kmer_vec[i];
			cnt = cnt_vec[i];

			Kmer tw = km.twin();
			Kmer rep = (km < tw) ? km : tw;
			
			if (kmap_ptr->contains(rep)) {
				tmp = kmap_ptr->find(rep) + cnt;
				if (tmp > 65535) tmp = 65535;
				kmap_ptr->update(rep, tmp);
			}
			//it = kmap_ptr->find(rep);
			//if (it != kmap_ptr->end()) {
			//	it->SetVal(it->GetVal()+cnt);
			//}
		}
	}
}

//record kmer counts in hash table for a single sample(kmc library)
float_t KmcCount(string kmc_file_name, ckhmap_t *kmap_ptr, bool verbose,  unsigned int threads_to_use, size_t ksize) {
//float_t KmcCount(string kmc_file_name, ckhmap_t *kmap_ptr, bool verbose, pool &tp, unsigned int num_threads) {
	CKMCFile kmer_data_base;

	if (!kmer_data_base.OpenForListing(kmc_file_name))
		exit(EXIT_FAILURE);
	else {
		CKMCFileInfo info;

		kmer_data_base.Info(info);

		//float counter;
		char str[ksize];
		uint32 cnt;
		float_t tot_coverage = 0;

		vector<Kmer> kmer_vec(info.total_kmers);
		vector<uint32> cnt_vec(info.total_kmers);
		
		//TODO: add min_count_to_set and max_count_to_set
		CKmerAPI kmer_object(info.kmer_length);

		uint64 idx_k = 0;
		while (kmer_data_base.ReadNextKmer(kmer_object, cnt)) {
			kmer_object.to_string(str);
			Kmer km(str);
			
			kmer_vec[idx_k] = km;
			cnt_vec[idx_k] = cnt;
			tot_coverage += log(cnt);
			idx_k ++;
		}

		/*pool tp(num_threads);
		if (verbose) {
			cout << "start " << num_threads << " threads" << endl;
		}*/
		#pragma omp parallel for num_threads(threads_to_use)
		for (unsigned int tid = 0; tid < threads_to_use; tid++) {
			SearchKmer(kmap_ptr, kmer_vec, cnt_vec, threads_to_use, tid);
			//tp.schedule(boost::bind(SearchKmer, kmap_ptr, boost::ref(kmer_vec), boost::ref(cnt_vec), num_threads, tid));
		}
		//tp.wait();
		//tp.clear();

		kmer_data_base.Close();
		kmer_vec.clear();
		cnt_vec.clear();

		return tot_coverage;
	
	}
}
