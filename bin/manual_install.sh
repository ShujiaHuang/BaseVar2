#!/bin/bash

# In MacOS
g++ -O3 -fPIC ../src/*.cpp ../src/io/*.cpp ../htslib/libhts.a -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -o basevar

# In Linux
g++ -O3 -fPIC ../src/*.cpp ../src/io/*.cpp ../htslib/libhts.a -I ../htslib -lz -lbz2 -lm -llzma -lpthread -lcurl -lssl -lcrypto -o basevar


