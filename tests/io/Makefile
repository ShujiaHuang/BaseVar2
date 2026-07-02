CXX = g++

LIBS = -lz -lm -lbz2 -llzma -lpthread -lcurl


HTSSRC = ../../htslib  # htslib
# Adjust $(HTSSRC) to point to your top-level htslib directory
ifdef HTSSRC
CPPFLAGS += -I"$(realpath $(HTSSRC))"
LIBS := "$(realpath $(HTSSRC))/libhts.a" $(LIBS)
else
$(info HTSSRC not defined, assuming systemwide installation)
LIBS += -lhts
endif

BSSRC = ../../src   # include all the header files(*.h) of BaseVar
CPPFLAGS += -I"$(realpath $(BSSRC))"

#modied from htslib makefile
FLAGS = -O3 -fPIC -std=c++17
FLAGS2 = $(CPPFLAGS) $(FLAGS) $(LDFLAGS)
CXXFLAGS := $(FLAGS2) $(CXXFLAGS)


#for compiling with ZSTD which is used for .bgen file format
ifeq ($(WITH_ZSTD),1)
LIBS += -lzstd
CXXFLAGS += -D__ZSTD__
endif

BSSRC = ../../src
CXXBS = $(wildcard $(BSSRC)/*.cpp)
CXXSRC = $(wildcard *.cpp)

OBJ = $(CXXSRC:.cpp=.o) $(CXXBS:.cpp=.o)


PROGRAMS = test
#all: $(PROGRAMS)


%.o: %.cpp    
	$(CXX) $(CXXFLAGS) -o $@ -c $<


test_fasta: test_fasta.o
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)
	#$(CXX) $(FLAGS) -o test_fasta  $(LIBS)

# `basevar motif` subcommand unit test.
# Links the whole BaseVar source tree (mirrors how the runner is built)
# minus main.cpp, which would conflict with the test's main().
MOTIF_BS_OBJ = $(BSSRC)/motif_counter.o $(BSSRC)/io/bam.o $(BSSRC)/io/bam_header.o \
               $(BSSRC)/io/bam_record.o $(BSSRC)/io/fasta.o $(BSSRC)/io/utils.o
test_motif_counter: test_motif_counter.o $(MOTIF_BS_OBJ)
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)

# `basevar fetalfrac` (hidden) subcommand unit test.
# Same linking strategy as test_motif_counter: pull in the module's
# implementation plus the BaseVar I/O abstractions it depends on.
FF_BS_OBJ = $(BSSRC)/fetal_fraction.o $(BSSRC)/io/bam.o $(BSSRC)/io/bam_header.o \
            $(BSSRC)/io/bam_record.o $(BSSRC)/io/utils.o
test_fetal_fraction: test_fetal_fraction.o $(FF_BS_OBJ)
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)

# `basevar caller` Bayesian genotype calling unit test.
# Phase 1: links algorithm.o (core Bayesian functions).
# Phase 2: adds caller_utils.o + basetype.o for pl_index_to_genotype and
#          process_sample_variant integration tests.
CALCPL_BS_OBJ = $(BSSRC)/algorithm.o $(BSSRC)/caller_utils.o $(BSSRC)/basetype.o \
                $(BSSRC)/io/utils.o $(BSSRC)/io/iobgzf.o $(BSSRC)/io/fasta.o
test_bayesian_caller: test_bayesian_caller.o $(CALCPL_BS_OBJ)
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)

# Binary batchfile (.bbf/.bbi) format unit test.
# Tests write/read roundtrip, index correctness, edge cases (long reads,
# long indels, empty positions, multi-sample, mixed variant types).
BBF_BS_OBJ = $(BSSRC)/io/batchfile_binary.o $(BSSRC)/io/iobgzf.o \
             $(BSSRC)/caller_utils.o $(BSSRC)/basetype.o \
             $(BSSRC)/algorithm.o $(BSSRC)/io/utils.o \
             $(BSSRC)/io/fasta.o
test_batchfile_binary: test_batchfile_binary.o $(BBF_BS_OBJ)
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)

# BGZF block-level utilities unit test.
# Header-only (iobgzf.h), only needs iobgzf.o for BGZFile class linkage.
BGZF_BS_OBJ = $(BSSRC)/io/iobgzf.o $(BSSRC)/io/utils.o
test_bgzf_utils: test_bgzf_utils.o $(BGZF_BS_OBJ)
	$(CXX) $(FLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f $(OBJ) test_fasta test_motif_counter test_motif_counter.tsv \
	      test_fetal_fraction test_fetal_fraction.tsv test_fetal_fraction.bed \
	      test_bayesian_caller test_batchfile_binary test_bgzf_utils


