#ifndef __HASH_LSH__
#define __HASH_LSH__

#include <algorithm>
#include <random>
#include <string>
#include <vector>



using namespace std;

namespace Hash{
typedef vector<float> hashFunction;
typedef vector<hashFunction> hashTable;
class LSH {
 public:
  static hashFunction generateNormalHashFunc(int size);
  static hashFunction generateUniformHashFunc(int size, float min, float max);
  static hashTable generateHashTable(int htable_size, int hfunction_size);

  // Different LSHs.
  static int random_projection(const vector<float> geneAbs, const hashFunction& function);
  static int random_projection(const vector<float> geneAbs, const hashTable& table);
  static string p_stable(const vector<float> geneAbs, const hashFunction& function, float b, float r);
  static string p_stable(const vector<float> geneAbs, const hashTable& table, float b, float r);
};

}  // namespace Core
#endif
