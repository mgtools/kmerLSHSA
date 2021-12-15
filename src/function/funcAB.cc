#include "funcAB.h"

namespace Core {

void AB::SetAbundance(Abundance* abundance, const vector<uint64_t>& ids, const vector<float>& values) {
  (*abundance)._values = values;
  (*abundance)._ids = ids;
}


/*
void AB::randAbundance(Abundance* abundance, const Abundance& ab, float scale ){
  random_device rd;
  mt19937 gen(rd());
  vector<float> new_gene_values;
  vector<float> gene_values = ab._values;
  vector<int> gene_locs = ab._locs;
  string title = ab._titles;
  int size = gene_values.size();
  int count = ab._count;
  float rand = 0.;
  normal_distribution<> dis(0,1);
  for(int i =0; i<size; ++i){
	rand = dis(gen) * scale;
	new_gene_values.push_back(rand + gene_values[i]);

  }
  SetAbundance(abundance, false, count, title, new_gene_values, gene_locs);
}
void AB::randomAbundance(Abundance* abundance, int count, float scale){
  vector<float> gene_ab;
  vector<int> gene_loc;
  float rand;
  string title = "default";
  random_device rd;
  mt19937 gen(rd());

  normal_distribution<> dis(0,1);
  for (int i=0; i<count; ++i){
	rand = dis(gen) * scale;
	if (rand > 0){
	  gene_ab.push_back(rand);
	  gene_loc.push_back(i);
	}
  }
  SetAbundance(abundance, false, false, count,1, title, title, gene_ab, gene_loc);
}
*/
void AB::SetConsensus(Abundance* abundance, const Abundance& ab1, const Abundance& ab2) {
  vector<float> new_ab_value;
  vector<float> ab_value1 = ab1._values, ab_value2 = ab2._values;
  vector<uint64_t> ab1_ids = ab1._ids, ab2_ids = ab2._ids;
  int ab1_count = ab1_ids.size(), ab2_count = ab2_ids.size();
  int all_count = ab1_count + ab2_count ;
  ab1_ids.insert(ab1_ids.end(),  ab2_ids.begin(), ab2_ids.end());
  /*
  cout << "AB1" << endl;
  cout << ab1 << endl;
  cout << "AB2" << endl;
  cout << ab2 << endl;
  */
  int size = ab_value1.size();
  new_ab_value.reserve(size);
  for(int i =0; i<size; i++){
    new_ab_value.push_back(ab_value1[i]*ab1_count/all_count+ab_value2[i]*ab2_count/all_count);
  }

  SetAbundance(abundance, ab1_ids, new_ab_value);
  //cout << "merged " << endl;
  //cout << (*abundance) << endl;
}

void AB::WRS(unordered_set<uint64_t> *group1, unordered_set<uint64_t> *group2, Abundance* abundance, int num_sample1, int num_sample2, float pvalue_thresh, int size_thresh){
  double freq_array1[num_sample1], freq_array2[num_sample2];
  double bothtails, lefttail, righttail;
  alglib::real_1d_array array_s1, array_s2;
  bothtails = 0;
  lefttail=0;
  righttail=0;
  unordered_set<uint64_t> kmers1, kmers2;
  swap(kmers1, *group1);
  swap(kmers2, *group2);
  vector<float> values = abundance->_values;
  vector<uint64_t> ids = abundance->_ids;
  int k = 0;
  if (ids.size() > size_thresh){
    for(int i=0; i<num_sample1; i++){
      freq_array1[i] = (double)values[i];
    }
    for(int j=0; j< num_sample2; j++){
      freq_array2[j] = (double)values[num_sample1+j];
    }

    array_s1.setcontent(num_sample1, freq_array1);
	  array_s2.setcontent(num_sample2, freq_array2);

	  //alglib::mannwhitneyutest(array_s1, num_sample1, array_s2, num_sample2, bothtails, lefttail, righttail);
    //alglib::ftest(array_s1, num_sample1, array_s2, num_sample2, bothtails, lefttail, righttail);
    alglib::studentttest2(array_s1, num_sample1, array_s2, num_sample2, bothtails, lefttail, righttail);

    if(lefttail <= pvalue_thresh) {
	    kmers2.insert(ids.begin(), ids.end());
	  } else if (righttail <= pvalue_thresh) {
	    kmers1.insert(ids.begin(), ids.end());
	  }
  }
  swap(*group1, kmers1);
  swap(*group2, kmers2);
}

void AB::dynamicST(unordered_set<uint64_t> *group1, unordered_set<uint64_t> *group2, vector<Abundance*>* clusteredab_ptr, int num_sample1, int num_sample2, float pvalue_thresh, int size_thresh, uint64_t extracted_size_min, bool verbose){
  float pval = pvalue_thresh;
  float pval_step = 0.0001;
  double freq_array1[num_sample1], freq_array2[num_sample2];
  double bothtails, lefttail, righttail;
  uint64_t extracted_size_thresh = extracted_size_min;
  alglib::real_1d_array array_s1, array_s2;
  bothtails = 0;
  lefttail=0;
  righttail=0;
  unordered_set<uint64_t> kmers1, kmers2;
  vector<Abundance*> unselected_ab, selected_ab1, selected_ab2, total_ab;
  swap(kmers1, *group1);
  swap(kmers2, *group2);
  swap(total_ab, *clusteredab_ptr);
  
  while((kmers1.size() <= extracted_size_thresh || kmers2.size() <= extracted_size_thresh) && pval <= 0.3){
    vector<Abundance*>().swap(unselected_ab);
    vector<Abundance*>().swap(selected_ab1);
    vector<Abundance*>().swap(selected_ab2);

    for(size_t i =0; i< total_ab.size(); i++){
      Abundance* abundance = total_ab[i];
      vector<float> values = abundance->_values;
      vector<uint64_t> ids = abundance->_ids;
      int k = 0;
      if (ids.size() > size_thresh){
        for(int i=0; i<num_sample1; i++){
          freq_array1[i] = (double)values[i];
        }
        for(int j=0; j< num_sample2; j++){
          freq_array2[j] = (double)values[num_sample1+j];
        }

        array_s1.setcontent(num_sample1, freq_array1);
	      array_s2.setcontent(num_sample2, freq_array2);

	      //alglib::mannwhitneyutest(array_s1, num_sample1, array_s2, num_sample2, bothtails, lefttail, righttail);
        //alglib::ftest(array_s1, num_sample1, array_s2, num_sample2, bothtails, lefttail, righttail);
        alglib::studentttest2(array_s1, num_sample1, array_s2, num_sample2, bothtails, lefttail, righttail);

        if(lefttail <= pval) {
          selected_ab2.push_back(abundance);
	      } 
        else if (righttail <= pval) {
          selected_ab1.push_back(abundance);
	      }
        else{
          unselected_ab.push_back(abundance);
        }
      }
    }

/*
    if(){

    }
    else if(){
      extracted_size_thresh = ;
    }
*/
    if(verbose){
      cout << "p-value threshold : " << pval << endl;
      cout << "The number of clusters that significantly enriched in Group A : " << selected_ab1.size() << endl;
      cout << "The number of clusters that significantly enriched in Group B : " << selected_ab2.size() << endl;
    }

    if (kmers1.size() <= extracted_size_thresh){
      for(size_t i =0; i< selected_ab1.size(); i++){
        vector<uint64_t> ids = selected_ab1[i]->_ids;
	      kmers1.insert(ids.begin(), ids.end());
      }
    }
    for(size_t i =0; i< selected_ab1.size(); i++){
      delete selected_ab1[i];
    }
    if (kmers2.size() <= extracted_size_thresh){
      for(size_t i =0; i< selected_ab2.size(); i++){
        vector<uint64_t> ids = selected_ab2[i]->_ids;
	      kmers2.insert(ids.begin(), ids.end());
      }
    }
    for(size_t i =0; i< selected_ab2.size(); i++){
      delete selected_ab2[i];
    }

    if(verbose){
      cout << "The number of k-mers enriched in Group A : " << kmers1.size() << endl; 
      cout << "The number of k-mers enriched in Group B : " << kmers2.size() << endl; 
    }

    swap(unselected_ab, total_ab);
    pval += pval_step;
  }
  swap(*clusteredab_ptr, total_ab);
  swap(*group1, kmers1);
  swap(*group2, kmers2);
}

void AB::SetMean(Abundance* abundance,const vector<Abundance*>& local_ab){
  int length = local_ab.size();
  Abundance* tmp = new Abundance();
  auto& previous = *(local_ab[0]);
  for(int i =1; i < length; ++i){
    //cout << "SetMean : " << i << endl;
    auto& current = *(local_ab[i]);
    SetConsensus(tmp, current, previous);
    previous = *tmp;
  }
  abundance = tmp;
}

bool AB::isSameAb(Abundance* ab1, Abundance* ab2){
  bool out = false;
  vector<float> ab_values1 = ab1->_values;
  vector<float> ab_values2 = ab2->_values;
  cout << "isSameAb start" << endl;
  
  for(int i=0; i<ab_values1.size(); ++i){
    if ( ab_values1[i] != ab_values2[i]){
      return out;
    }
  }
  out = true;
  return out;
}
}
