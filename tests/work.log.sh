
g++ -O3 -fPIC ../src/main.cpp ../src/basetype.h ../src/basetype.cpp ../src/basetype_caller.cpp ../src/utils.cpp ../src/fasta.cpp ../src/bam_header.cpp ../src/bam.cpp ../src/bam_record.cpp ../htslib/libhts.a -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o basevar


