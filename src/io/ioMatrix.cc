#include "ioMatrix.h"

namespace Utility {

int IOMat::parseLine(char* line){
  //This assumes that a digit will be found and the line ends in " Kb".
  int i = strlen(line);
  const char* p = line;
  while (*p <'0' || *p > '9') p++;
  line[i-3] = '\0';
  i = atoi(p);
  return i;
}

int IOMat::getValue(int i){ //Note: this value is in KB!
  FILE* file = fopen("/proc/self/status", "r");
  int result = -1;
  char line[128];

  while (fgets(line, 128, file) != NULL){
    if (strncmp(line, "VmSize:", 7) == 0){
      result = parseLine(line);
      break;
    }
  }
  fclose(file);
  cout << i << ". Virtual Memory Usage(KB) : " << result << endl;
  return result;
}

float IOMat::convert(char const* source, char ** endPtr ) {
  char* end;
  int left = strtol( source, &end, 10 );
  float results = left;
  if ( *end == '.' ) {
      char* start = end + 1;
      int right = strtol( start, &end, 10 );
      static float const fracMult[]
          = { 0.0, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001 };
      results += right * fracMult[ end - start ];
  }
  if ( endPtr != nullptr ) {
      *endPtr = end;
  }
  return results;
}

void IOMat::ReadClusterAll(vector<Abundance*>* AbundanceMat, int num_samples, string file_name, bool verbose){
  auto start_time = chrono::high_resolution_clock::now();
  string clust_file_name = file_name +".clust" ;
  uint64_t loc = 0;
  uint64_t id;
  string line;
  const char* lineStart;
  char* lineEnd;

  size_t line_cnt;
  
  /////read binary file
  ifstream inbfile(file_name, ios::in | ios::binary);
  if (!inbfile.is_open()) {
    cout << "Error! file ( " << file_name << " ) not open!" << endl;
    exit(-1);
  }
  
  inbfile.seekg(0, ios::end);
  line_cnt = inbfile.tellg() / (sizeof(float)*num_samples);
  inbfile.seekg(0, ios::beg);
  cout << "read start : "<< file_name << ", size :" << line_cnt << endl;
  
  vector<vector<float>> values(line_cnt, vector<float>(num_samples));
  vector<Abundance*> local_abundance(line_cnt);

  for(size_t i =0; i < line_cnt; i++){
    //read gene Name
    inbfile.read(reinterpret_cast<char*> (values[i].data()), sizeof(float)*num_samples);
  }
  inbfile.close();

  //////read cluster file
  ifstream infile(clust_file_name);
  if (!infile.is_open()) {
    cout << "Error! file ( " << clust_file_name << " ) not open!" << endl;
    exit(-1);
  }
  cout << "read cluster result : "<< clust_file_name << endl;
  
  while(getline(infile, line)){
    //ids.clear();

    lineStart = line.c_str();
    id = strtol(lineStart, &lineEnd, 10);
	  vector<uint64_t> ids(id);
	  uint64_t numID = 0;
    while(lineStart != lineEnd && numID < ids.size()){
      // https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-floats-in-c-quickly
      //ids.push_back(id);
      lineStart = lineEnd;
      id =strtol(lineStart, &lineEnd, 10);
	    ids[numID] = id;
	    numID++;
    }
    //AB::UpdateAbundanceIDs(AbundanceMat->at(loc), new_ids);
    Abundance* abundance = new Abundance();
    AB::SetAbundance(abundance, ids, values[loc]);
    local_abundance[loc] = abundance;
	  loc++;
	  vector<uint64_t>().swap(ids);
    infile.clear();
  }
  infile.close();
  swap(*AbundanceMat, local_abundance);
  vector<vector<float>>().swap(values);
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time).count();
  if(verbose){
    cout << "Reading the result file takes secs: " << elapsed_read << endl;
  }
}

void IOMat::ReadCluster(vector<Abundance*>* AbundanceMat, int num_samples, string file_name,  streamoff start_line, uint64_t num_lines, bool verbose){
  auto start_time = chrono::high_resolution_clock::now();
  string clust_file_name = file_name +".clust" ;
  uint64_t loc = 0;
  uint64_t id;
  string line;
  const char* lineStart;
  char* lineEnd;
  //vector<int> ids;
  
  /////read binary file
  ifstream inbfile(file_name, ios::in | ios::binary);
  if (!inbfile.is_open()) {
    cout << "Error! file ( " << file_name << " ) not open!" << endl;
    exit(-1);
  }
  inbfile.seekg( start_line * sizeof(float) * num_samples , ios::beg);
  cout << "read start : "<< file_name << ", size :" << num_lines << endl;
  
  vector<vector<float>> values(num_lines, vector<float>(num_samples));
  vector<Abundance*> local_abundance(num_lines);

  for(size_t i =0; i < num_lines; i++){
    //read gene Name
    inbfile.read(reinterpret_cast<char*> (values[i].data()), sizeof(float)*num_samples);
  }
  inbfile.close();

  //////read cluster file
  ifstream infile(clust_file_name);
  if (!infile.is_open()) {
    cout << "Error! file ( " << clust_file_name << " ) not open!" << endl;
    exit(-1);
  }
  cout << "read cluster result : "<< clust_file_name << endl;
  
  for(size_t i =0; i< start_line; i++){
    //getline(infile, line);
	  infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  while(getline(infile, line) && loc < num_lines){

    lineStart = line.c_str();
    id = strtol(lineStart, &lineEnd, 10);
	  //cout << "size : " << id << endl;
	  vector<uint64_t> ids(id);
	  uint64_t numID = 0;

    while(lineStart != lineEnd && numID < ids.size()){
      // https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-floats-in-c-quickly
      //ids.push_back(id);
      lineStart = lineEnd;
      id =strtol(lineStart, &lineEnd, 10);
	  //cout << "numID : " << numID << " id : " << id  <<endl;
	    ids[numID] = id;
	    numID++;
      infile.clear();
    }
    //AB::UpdateAbundanceIDs(AbundanceMat->at(loc), new_ids);
    Abundance* abundance = new Abundance();
    AB::SetAbundance(abundance, ids, values[loc]);
    local_abundance[loc] = abundance;
	  //cout << *abundance << endl;
	  loc++;
	  vector<uint64_t>().swap(ids);
  }
  infile.close();
  swap(*AbundanceMat, local_abundance);
  vector<vector<float>>().swap(values);
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time).count();
  if(verbose){
    cout << "Reading the result file takes secs: " << elapsed_read << endl;
  }
}

void IOMat::ReadMatrix(vector<Abundance*>* AbundanceMat, bool normalization,  string file_name ) {
  int sample_cnt = 0;
  uint64_t line_cnt = 0;
  float totalvalueAb= 0.0;
  float valueAb;
  vector<float> nonzero_ab;
  vector<int> nonzero_loc;
  vector<uint64_t> ids;
  string line, head;
  ifstream infile(file_name);
  if (!infile.is_open()) {
    cout << "Error! file ( " << file_name << " ) not open!" << endl;
    exit(-1);
  }

  cout << "read start : "<< file_name << endl;
  const char* lineStart;
  char* lineEnd;

  while (getline(infile, line)){
    nonzero_ab.clear();
    nonzero_loc.clear();
    ids.clear();
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (line[0] == '\t'){
      head = line;
      continue;
    }
    //read gene Name

    lineStart = line.c_str();
    valueAb = strtod(lineStart, &lineEnd);
    ids.push_back(line_cnt);

    while(lineStart != lineEnd){
      // https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-floats-in-c-quickly
      
      nonzero_ab.push_back(valueAb);
      if(normalization){totalvalueAb += pow(valueAb, 2);}
      
	    ++sample_cnt;
	    lineStart = lineEnd;
      valueAb = strtod(lineStart, &lineEnd);
    }

    //check using normalization
    if (normalization){
      for(int i = 0; i< nonzero_ab.size(); i++){
        nonzero_ab[i] /= sqrt(totalvalueAb);
      }
    }

    Abundance* abundance = new Abundance();
    AB::SetAbundance(abundance, ids, nonzero_ab);
	//cout << *abundance << endl;
    (*AbundanceMat).push_back(abundance);
    ++line_cnt;
    sample_cnt = 0;
    totalvalueAb = 0;
  }
  cout << "read #total lines: " << line_cnt << endl;
  infile.close();
}


void IOMat::SaveResult(const vector<Abundance*>* abs_all, string out_file, bool delfile, int ignore_small, bool verbose){
  vector<uint64_t> ids;
  auto start_time = chrono::high_resolution_clock::now();
  if (delfile){
    if(remove(out_file.c_str()) != 0){
      perror("File deletion failed");
    }
    else{
      cout << "files are removed" << endl;
    }
  }
  ofstream out;
  out.open(out_file, fstream::out | fstream::app);
  for (size_t i = 0; i < abs_all->size(); ++i) {
	  ids = abs_all->at(i)->_ids;
    if(ids.size()> ignore_small){
      out << ids.size();
		  for(int j=0; j< ids.size(); j++){
			  out << "\t" << ids[j] ;
		  }
		  out<< endl;
	  }
  }
  out.close();
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time).count();
  if(verbose){
    cout << "Saveing the result file takes secs: " << elapsed_read << endl;
  }
}

void IOMat::SaveMatrix(const vector<Abundance*>* abs_all, string out_file, bool delfile, int ignore_small) {
  if (delfile){
    if(remove(out_file.c_str()) != 0){
      perror("File deletion failed");
    }
    else{
      cout << "files are removed" << endl;
    }
  }
  ofstream out;
  out.open(out_file, fstream::out | fstream::app);

	for (size_t i = 0; i < abs_all->size(); ++i) {
		vector<float> values = abs_all->at(i)->_values;
    size_t idsSize = (abs_all->at(i)->_ids).size();
    int dim = values.size();
    if(idsSize > ignore_small){
		  for(int j=0; j<dim-1; ++j){
			  out << values[j] << "\t";
      }	
		  out << values[dim-1] << endl;
    }
	}
	out.close();
}

void IOMat::SaveBinary(const vector<Abundance*>* abs_all, string out_file, bool delfile, int ignore_small, bool verbose) {
  auto start_time = chrono::high_resolution_clock::now();
  if (delfile){
    if(remove(out_file.c_str()) != 0){
      perror("File deletion failed");
    }
    else{
      cout << "files are removed" << endl;
    }
  }
  ofstream outfile;
  outfile.open(out_file, ios::out | ios::binary | ios::app);

	for (size_t i = 0; i < abs_all->size(); ++i) {
		vector<float> values = abs_all->at(i)->_values;
    int idsSize = (abs_all->at(i)->_ids).size();
    int dim = values.size();
    if(idsSize > ignore_small){
		  //for(int j=0; j<dim; ++j){
			outfile.write(reinterpret_cast<const char*>( &values[0] ), sizeof( float )*dim) ;
      //}	
    }
	}
	outfile.close();
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time).count();
  if(verbose){
    cout << "Saveing the binary file takes secs: " << elapsed_read << endl;
  }
}

void IOMat::convertHTMat(uint16_t **ary_count, vector<float_t> &v_kmers, int tot_sample, bool verbose, uint64_t batch_size, streamoff batch_offset, vector<Abundance*>* unknown_abundance_ptr){
  vector<Abundance*> abVec, *abVec_ptr;
  abVec_ptr = &abVec;
  uint64_t cnt, total_cnt;
  float_t num_kmers;
  vector<float> values;
  vector<uint64_t> ids;
  float value, totalValue;
  auto start_time = chrono::high_resolution_clock::now();

  totalValue = 0.0;
  total_cnt = 0;

  ids.reserve(1);
  values.reserve(tot_sample);
  abVec.reserve(batch_size);

  //swap(abVec, *unknown_abundance_ptr);

  for(uint64_t i=0; i<batch_size; i++){
    ids.push_back(batch_offset+i);
    for (int j =0; j<tot_sample;j++){
      cnt = ary_count[j][i];
      num_kmers = v_kmers[j];
	  total_cnt += cnt;
      value = float(log(cnt+1.0)) - num_kmers ;
      values.push_back(value);
    }
    if(total_cnt >0.01* tot_sample){
      Abundance *abundance = new Abundance();
      AB::SetAbundance(abundance, ids, values);
  	//cout << *abundance << endl;
	  abVec.push_back(abundance);
	 }

    totalValue = 0.0;
	total_cnt = 0;
    values.clear();
    ids.clear();
  }
  //cout << "print : " << abVec.size() << endl;
  //unknown_abundance_ptr->insert(unknown_abundance_ptr->end(), abVec.begin(), abVec.end());
  swap(*unknown_abundance_ptr, abVec);
  //SaveMatrix(abVec_ptr, file_name, tot_sample, false);
  //SaveResult(abVec_ptr, file_name+".clust", false);
  for(int k=0; k<abVec.size(); k++){
    delete abVec[k];
  }
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<float>>(end_time - start_time).count();
  
  if(verbose){
    cout << "Finish load partial matrix" << endl;
	  cout << "loading Matrix takes secs:\t" << elapsed_read << endl;
  }
}

}  // namespace Utility
