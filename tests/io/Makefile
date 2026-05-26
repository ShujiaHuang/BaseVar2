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

clean:
	rm -f $(OBJ) test_fasta test_motif_counter test_motif_counter.tsv


