#!/bin/sh

g++ -c -fopenmp -lgomp -fPIC rngstream.cpp

PKG_CXXFLAGS="-c -fopenmp -w -std=c++11"  R CMD SHLIB SMMB_revised.cpp -fopenmp rngstream.o

PKG_CXXFLAGS="-c -fopenmp -w -std=c++11"  R CMD SHLIB SMMB_source.cpp -fopenmp rngstream.o