g++ -O3 -fPIC test_fasta.cpp ../../src/fasta.cpp ../../src/utils.cpp ../../htslib/libhts.a -I ../../src -I ../../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_fasta && ./test_fasta

make
