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
FLAGS = -O3 -fPIC
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


PROGRAMS = test_fasta

all: $(PROGRAMS)


%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $*.cpp


$(PROGRAMS): $(OBJ)
	$(CXX) $(FLAGS) -o test_fasta *.o $(LIBS)

clean:
	rm -f *.o


