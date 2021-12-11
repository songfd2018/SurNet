#!/bin/sh

g++ -c -fopenmp -lgomp -fPIC rngstream.cpp

PKG_CXXFLAGS="-c -fopenmp -w -std=c++11"  R CMD SHLIB SMMB_cpp.cpp -fopenmp rngstream.o
