ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit -lFoam

CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./
LIBS           = $(ROOTLIBS) 

GLIBS          = $(filter-out -lNew, $(NGLIBS))

CXX	      += -I./lib/
OUTLIB	      = ./lib/
.SUFFIXES: .cc,.C
.PREFIXES: ./lib/

#----------------------------------------------------#

main: lib/digit_b2.o

lib/digit_b2.o:
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)digit_b2.o $<
lib/digit_b2.o: src/digit_b2.cc

clean:
	rm -f src/*~
	rm -f *~
	rm -f $(OUTLIB)*.o

#----------------------------------------------------#
