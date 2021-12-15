#include "cluster.h"


void p_lsh(vector<Abundance*> *hash_values, vector<int> *hash_keys, hashTable& hash_table, Abundance* abundance) {
  
  //vector<Abundance*> local_hash_values;
  //vector<int> local_hash_keys;
  int key = LSH::random_projection(abundance->_values, hash_table);
  hash_values->push_back(abundance);
  hash_keys->push_back(key);
  
}


void merge_hashtable(vector<vector<Abundance*>>* final_table, size_t hash_func_num, const vector<vector<Abundance*>*>& part_hash_values,  const vector<vector<int>*>& part_hash_keys) {
  const size_t buckets = size_t(pow(2, hash_func_num));
  vector<vector<Abundance*>> local_table(buckets, vector<Abundance*>());

  for (int i = 0; i < part_hash_keys.size(); ++i) {
    for (int j = 0; j < part_hash_keys[i]->size(); ++j) {
      const auto& key = part_hash_keys[i]->at(j);
      const auto& ptr = part_hash_values[i]-> at(j);
      local_table[key].push_back(ptr);
    }
  }
  final_table->swap(local_table);
}

void free_hashtable(vector<vector<Abundance*>>* final_table, int s){
  size_t size = (final_table->at(s)).size();
  for(size_t i = 0; i< size; ++i){
	  (final_table->at(s)).pop_back();
  }
  vector<Abundance*>().swap(final_table->at(s));
}
void merge_abundance(vector<Abundance*>* final_abundance, const vector<vector<Abundance*>*>& part_abundance) {
  vector<Abundance*> local_abundance;
  for (const auto& abundance : part_abundance) {
    local_abundance.insert(local_abundance.end(), abundance->begin(), abundance->end());
  }
  final_abundance->swap(local_abundance);
}

void free_abundance(vector<Abundance*> final_abundance){
  vector<Abundance*> local_abundance;
  for (size_t i=0; i< final_abundance.size(); ++i){
	  delete final_abundance[i];
  }
  final_abundance.clear();
  vector<Abundance*>().swap(final_abundance);
}

void p_cluster(vector<Abundance*> *part_ab, vector<Abundance*>* candidates, float threshold) {
  float distance;
  size_t size;
  vector<Abundance*> local_ab;
  size = candidates->size();
  size_t i = 1;
  size_t j = 0;
  while(i < size) {
    auto& current = *(candidates->at(i));
    
    for (j = 0; j < i; ++j) {
      auto& candidate = *(candidates->at(j));
      distance = Distance::cosine(current._values , candidate._values);
      if (1 - distance >= threshold) {
        Abundance* ab_new = new Abundance();
        AB::SetConsensus(ab_new, current, candidate);
        delete candidates->at(i);
        delete candidates->at(j);
        candidates->at(j) = ab_new;
        candidates->at(i) = candidates->at(--size);
        candidates->at(size) = NULL;
        break;
      }
    }
    if (j == i){
      ++i;
    }
  }
  part_ab->insert(part_ab->end(), candidates->begin(), (candidates->begin())+size);
  
  //swap(*part_ab, local_ab);
}

void nestedCluster(vector<Abundance*>* unknown_abundance_ptr, float similarity, int dim, int threads_to_use, bool verbose) {
  size_t unknown_abundance_size = 0;

  
  auto start_time_total = chrono::high_resolution_clock::now(),
  start_time = chrono::high_resolution_clock::now(),
  end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time).count();

  size_t hash_func_num = floor(log2((*unknown_abundance_ptr).size()));
  size_t buckets = size_t(pow(2, hash_func_num));
  vector<vector<Abundance*>> lsh_table(buckets, vector<Abundance*>());
  int abundance_size = 0;

  abundance_size = (*unknown_abundance_ptr).size();
    
  //if (verbose){
  //  cout << "1-iter clustering cos sim threshold:\t" << similarity <<" dimension : " <<dim<< endl;
  //  cout << "Size of profilings : " << (*unknown_abundance_ptr).size() << endl;
//	}
  
  hashTable hash_table = LSH::generateHashTable(hash_func_num, dim);

  // Apply LSH.
  auto start_time_iteration = chrono::high_resolution_clock::now();
  auto start_time_hashing = chrono::high_resolution_clock::now();

  vector<vector<Abundance*>*> part_hash_values;
  vector<vector<int>*> part_hash_keys;
  for (int i = 0; i < threads_to_use; ++i) {
    part_hash_values.push_back(new vector<Abundance*>());
    part_hash_keys.push_back(new vector<int>());
  }
    
  //pool tp(threads_to_use);
  #pragma omp parallel for num_threads(threads_to_use)
  for (size_t i = 0; i < abundance_size; ++i) {
    int tid = omp_get_thread_num();
    p_lsh( part_hash_values[tid],part_hash_keys[tid], hash_table, unknown_abundance_ptr->at(i));
  }
  merge_hashtable(&lsh_table, hash_func_num, part_hash_values, part_hash_keys);   
 
  for (int i = 0; i < threads_to_use; ++i) {
    delete part_hash_values[i];
    delete part_hash_keys[i];
  }

  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time_hashing).count();
  //if(verbose){
    //cout << "nested clustering hashing takes secs:\t" << elapsed_read << endl;
  //}

  auto start_time_clustering = chrono::high_resolution_clock::now();

  // Cluster within LSH buckets.
  
  vector<vector<Abundance*>*> part_abundance;
  for (int i =0; i< threads_to_use; ++i){
    part_abundance.push_back(new vector<Abundance*>());
  }
  
  int buckets_per_thread = lsh_table.size() / threads_to_use;
  
  #pragma omp parallel for num_threads(threads_to_use) schedule(dynamic)
  for (size_t s = 0; s < buckets; ++s) {
    auto& candidates = lsh_table[s];
    size_t size = candidates.size();
    int tid = omp_get_thread_num();
    p_cluster(part_abundance[tid], &candidates, similarity);
  }

  for (size_t i=0; i<buckets; i++){
    vector<Abundance*>().swap(lsh_table[i]);
  }
    
  auto start_time_merge = chrono::high_resolution_clock::now();
  merge_abundance(unknown_abundance_ptr, part_abundance);

  for(int i=0; i<threads_to_use; ++i){
    delete part_abundance[i];
  }
  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time_merge).count();
  //if(verbose){
  //  cout << "nested merging takes secs:\t" << elapsed_read << endl;
  //  cout << "nested #k-mers after clustering:\t" << (*unknown_abundance_ptr).size() << endl;
  //  cout << endl;
  //}
}

// get the lsh and iterations.
void Cluster(vector<Abundance*>* unknown_abundance_ptr, float min_similarity, float max_similarity, int cluster_iteration, unsigned int threads_to_use, int dim, int bucket_size_threshold, bool verbose) {
  int unknown_abundance_size = 0;
  string head = "";

  auto start_time_total = chrono::high_resolution_clock::now(),
  start_time = chrono::high_resolution_clock::now(),
  end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time).count();

  //float max_similarity = 0.95;  // Heuristics for maximum spearman similarity.
  float sim_step = (max_similarity - min_similarity)/cluster_iteration;
  float threshold = max_similarity;
  int iter = 0;
  size_t hash_func_num = floor(log2((*unknown_abundance_ptr).size()));
  size_t buckets = int(pow(2, hash_func_num));
  vector<vector<Abundance*>> lsh_table(buckets, vector<Abundance*>());
  size_t abundance_size = 0;

  while (iter++ < cluster_iteration) {
    abundance_size = (*unknown_abundance_ptr).size();
    
    if (iter > 1 ){
      hash_func_num = floor(log2((*unknown_abundance_ptr).size()));
      buckets = size_t(pow(2, hash_func_num));  
      lsh_table.resize(buckets);
    }
    
	  if (verbose){
      //cout << "Iteration:\t" << iter << ", cos sim threshold:\t" << threshold <<" dimension : " <<dim<< endl;
	  // cout << "Size of profilings : " << (*unknown_abundance_ptr).size() << endl;
	  }
    // Generate LSH.
    hashTable hash_table = LSH::generateHashTable(hash_func_num, dim);

    // Apply LSH.
    auto start_time_iteration = chrono::high_resolution_clock::now();
    auto start_time_hashing = chrono::high_resolution_clock::now();

    if(verbose){
      IOMat::getValue(1);
    }
    //vector<vector<Abundance*>> part_hash_values(threads_to_use, vector<Abundance*>());
    //vector<vector<int>> part_hash_keys(threads_to_use, vector<int>());
    vector<vector<Abundance*>*> part_hash_values;
    vector<vector<int>*> part_hash_keys;
    for (int i = 0; i < threads_to_use; ++i) {
      part_hash_values.push_back(new vector<Abundance*>());
      part_hash_keys.push_back(new vector<int>());
    }
    
    #pragma omp parallel for num_threads(threads_to_use)
    for (int i = 0; i < abundance_size; ++i) {
      int tid = omp_get_thread_num();
      p_lsh( part_hash_values[tid], part_hash_keys[tid], hash_table, unknown_abundance_ptr->at(i));
    }
    
    // Merge thread results.
    //if(verbose){
      //IOMat::getValue(2);
      //cout << "Merge thread results " << endl;
    //}

    merge_hashtable(&lsh_table, hash_func_num, part_hash_values, part_hash_keys);
    
    //if(verbose){    
      //IOMat::getValue(3);
    //}

    for (int i = 0; i < threads_to_use; ++i) {
      delete part_hash_values[i];
      delete part_hash_keys[i];
    }
    
    //swap(*unknown_abundance_ptr, vector<Abundance*> ());

    end_time = chrono::high_resolution_clock::now();
    elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time_hashing).count();
    //if(verbose){
      //IOMat::getValue(4);
      //cout << "hashing takes secs:\t" << elapsed_read << endl;
    //}

    auto start_time_clustering = chrono::high_resolution_clock::now();

    // Cluster within LSH buckets.
    //vector<vector<Abundance*>> part_abundance(threads_to_use, vector<Abundance*>());
    vector<vector<Abundance*>*> part_abundance;
    for (int i =0; i< threads_to_use; ++i){
      part_abundance.push_back(new vector<Abundance*>());
    }
  
    int buckets_per_thread = lsh_table.size() / threads_to_use;

    omp_set_nested(1);
    omp_set_max_active_levels(2);
    
    #pragma omp parallel for num_threads(threads_to_use) schedule(dynamic)
    for (size_t s = 0; s < buckets; ++s) {
      auto& candidates = lsh_table[s];
      size_t size = candidates.size();
      int tid = omp_get_thread_num();
      if(size > bucket_size_threshold){
        nestedCluster(&candidates, threshold, dim, 3, verbose);
        part_abundance[tid]->insert(part_abundance[tid]->end(), candidates.begin(), candidates.end());
      }
      else{
        p_cluster(part_abundance[tid], &candidates, threshold);
      }
    }
    omp_set_nested(0);

    for (int i=0; i<buckets; i++){
      vector<Abundance*>().swap(lsh_table[i]);
    }

    end_time = chrono::high_resolution_clock::now();
    elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time_clustering).count();
    if(verbose){
      //IOMat::getValue(6);
      cout << "clustering takes secs:\t" << elapsed_read << endl;
      cout << endl;
    }
    ////////////
    auto start_time_merge = chrono::high_resolution_clock::now();
    merge_abundance(unknown_abundance_ptr, part_abundance);

    for(int i=0; i<threads_to_use; ++i){
      delete part_abundance[i];
    }
    IOMat::getValue(8);

    end_time = chrono::high_resolution_clock::now();
    elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time_merge).count();
    if(verbose){
      cout << "merging takes secs:\t" << elapsed_read << endl;
      cout << "#k-mers after clustering:\t" << (*unknown_abundance_ptr).size() << endl;
      cout << endl;
    }

    threshold -= sim_step;
  }


  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time_total).count();
  if(verbose){
    cout << "kmerLSH algorithm hash+cluster takes (secs): " << elapsed_read << endl;
  }
  
}



