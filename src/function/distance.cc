#include "distance.h"
namespace Core{
/*
float Distance::spearsman(const Abundance& lhs, const Abundance& rhs){
   float distance = 0.
   alglib::real_1d_array array_s1, array_s2;
   arr
   array_s1.setcontent(num_sample1, freq_array1);
   array_s2.setcontent(num_sample2, freq_array2);
   alglib::mannwhitneyutest(array_s1, num_sample1, array_s2, num_sample2, bothtails, lefttail, righttail);

   //distance = 


}
*/
//TODO: Implement a highly efficient method like 'cosine'.

float Distance::euclidean(const vector<float> lhs,  const vector<float> rhs) {
  float distance = 0;
  for(int i=0; i < rhs.size(); i++){
    distance += pow( rhs[i] - lhs[i] ,2);
  }
  return sqrt(distance);
}

float Distance::cosine(const vector<float> lhs,  const vector<float> rhs) {
  float similarity = 0;
  float magnitude_lhs = 0, magnitude_rhs = 0;
  //cout << "cosine lhs : "<< lhs.size() << " rhs : " << rhs.size() << endl;
	for (int i = 0; i<lhs.size(); i++){
    similarity += lhs[i] * rhs[i];
    magnitude_lhs += lhs[i] * lhs[i]; 
    magnitude_rhs += rhs[i] * rhs[i];
	}
  similarity /= sqrt(magnitude_lhs) * sqrt(magnitude_rhs);
  return 1 - similarity;
}
  


}  //namespace Core
