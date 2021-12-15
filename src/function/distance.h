#ifndef __DISTANCE_H__
#define __DISTANCE_H__

#include <cmath>
#include <utility>
#include <unordered_map>
#include <vector>


#include <iostream>
using namespace std;


using std::pair;
using std::unordered_map;
using std::vector;

namespace Core {
class Distance {
 public:
  static float euclidean(const vector<float> lhs, const vector<float> rhs);
  static float cosine(const vector<float> lhs, const vector<float> rhs);
  //static float spearsman(const GeneAbs& lhs, const GeneAbs& rhs);

};

}  // namespace Core
#endif
