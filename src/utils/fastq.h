#ifndef FASTQ_H
#define FASTQ_H

#include <stdio.h>
#include <string.h>
#include <zlib.h>


#include <vector>
#include <string>

#include "../common/Common.h"
#include "../kmer/kseq.h"

// initialize kseq structures
//KSEQ_INIT(gzFile, gzread);

struct ReadEntry {
	char name[8196];
	char s[8196];
	char qual[8196];
	size_t name_len;
	size_t len;
};

class FastqFile {
 public:
  FastqFile(const vector< string> fnames);

  ~FastqFile();

  void close();
  void reopen();
  int read_next(char *read, size_t *read_len, char *seq, size_t *seq_len, unsigned int *file_id, char* qual = 0);
  int read(vector<ReadEntry> &read_vec, size_t &num_reads);
  uint64_t part_size = 1 << 16;


 private:
  vector<string>::const_iterator open_next();

  vector<string>::const_iterator fnit;
  unsigned int file_no;
  vector<string> fnames;
  gzFile fp;
  kseq_t *kseq;
};

#endif
