g++ -O3 -fPIC test_fasta.cpp ../../src/fasta.cpp ../../src/utils.cpp ../../htslib/libhts.a -I ../../src -I ../../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_fasta && ./test_fasta

g++ -O3 -fPIC test_bamheader.cpp ../../src/bam_header.cpp ../../src/utils.cpp ../../htslib/libhts.a -I ../../src -I ../../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_bamheader && ./test_bamheader

g++ -O3 -fPIC test_bamrecord.cpp ../../src/bam_record.cpp ../../src/bam_header.cpp ../../src/utils.cpp ../../htslib/libhts.a -I ../../src -I ../../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_bamrecord && ./test_bamrecord


g++ -O3 -fPIC test_bam.cpp ../../src/bam.cpp ../../src/bam_header.cpp ../../src/bam_record.cpp ../../src/utils.cpp ../../htslib/libhts.a -I ../../src -I ../../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o test_bam && ./test_bam

make
