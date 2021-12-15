
KMER_DIR := src/kmer
FUNC_DIR := src/function
UTIL_DIR := src/utils
HASH_DIR := src/hash
IO_DIR := src/io
KMC_API_DIR := $(KMER_DIR)/kmc_api
ALGLIB_DIR := $(UTIL_DIR)/alglib-3.17.0/src
CUR_DIR := $(notdir $(shell pwd))

INC := -I ./
LIBS := -L ./
CXX := g++
CXXFLAGS := -std=c++11 -O3 $(INC) -fopenmp 
CXXFLAGS1 := -std=c++11 -O3 $(INC)
CXXFLAGS2 := -std=c++11 -O3 
//CXXLINK := -O3 -std=c++11 -lm -lz -lpthread -lboost_system -lboost_thread
CXXLINK := -O3 -std=c++11 -lm -lz -lpthread -fopenmp

KMER_OBJS := $(KMER_DIR)/Kmer.o $(KMER_DIR)/KmerIntPair.o $(KMER_DIR)/kmc_reader.o
UTIL_OBJS := $(UTIL_DIR)/fastq.o
HASH_OBJS := $(HASH_DIR)/hash.o $(HASH_DIR)/lshash.o
IO_OBJS := $(IO_DIR)/ioHT.o $(IO_DIR)/ioMatrix.o $(IO_DIR)/ioFastQ.o
FUNC_OBJS := $(FUNC_DIR)/cluster.o $(FUNC_DIR)/distance.o $(FUNC_DIR)/funcAB.o

KMC_API_OBJS := \
$(KMC_API_DIR)/mmer.o \
$(KMC_API_DIR)/kmc_file.o \
$(KMC_API_DIR)/kmer_api.o

ALGLIB_OBJS := \
$(ALGLIB_DIR)/ap.o \
$(ALGLIB_DIR)/alglibinternal.o \
$(ALGLIB_DIR)/alglibmisc.o \
$(ALGLIB_DIR)/linalg.o \
$(ALGLIB_DIR)/specialfunctions.o \
$(ALGLIB_DIR)/statistics.o

$(KMER_DIR)/%.o: $(KMER_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@ 
$(UTIL_DIR)/%.o: $(UTIL_DIR)/%.cc
	$(CXX) $(CXXFLAGS1) -c $< -o $@
$(HASH_DIR)/%.o: $(HASH_DIR)/%.cc
	$(CXX) $(CXXFLAGS1) -c $< -o $@
$(IO_DIR)/%.o: $(IO_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(FUNC_DIR)/%.o: $(FUNC_DIR)/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@
$(KMC_API_DIR)/%.o: $(KMC_API_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ 
$(ALGLIB_DIR)/%.o: $(ALGLIB_DIR)/%.cpp
	$(CXX) $(CXXFLAGS2) -c $< -o $@

bin/kmerLSHSA: src/app/kmerLSH.cc $(KMER_OBJS) $(UTIL_OBJS) $(HASH_OBJS) $(IO_OBJS) $(FUNC_OBJS) $(KMC_API_OBJS) $(ALGLIB_OBJS)
	$(CXX) -o $@ $^ $(INC) $(LIBS) $(CXXLINK)

clean: 
	-rm -f $(KMER_DIR)/*.o
	-rm -f $(HASH_DIR)/*.o
	-rm -f $(IO_DIR)/*.o
	-rm -f $(UTIL_DIR)/*.o
	-rm -f $(FUNC_DIR)/*.o
	-rm -f $(KMC_API_DIR)/*.o
	-rm -f $(ALGLIB_DIR)/*.o
	-rm -f bin/kmerLSHSA

all: bin/kmerLSHSA




