################################################################################
#
# PhysicsTools/TagAndProbe Makefile (for standalone use outside CMSSW/SCRAM)
# ---------------------------------------------------------------------------
#
# adapted from SusyAnalysis/SusyNtuplizer/jec/JetMETObjects/Makefile
################################################################################

ifeq ($(STANDALONE_DIR),)
	standalone_dir:=..
	export STANDALONE_DIR:=$(standalone_dir)
endif


TMPDIR       = $(STANDALONE_DIR)/tmp
LIBDIR       = $(STANDALONE_DIR)/lib
BINDIR       = $(STANDALONE_DIR)/bin



CXX          = g++


ROOTCXXFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
INC          = $(shell $(ROOTSYS)/bin/root-config --incdir)
CXXFLAGS     = -O3 -Wall -fPIC -DSTANDALONE -I../.. -I/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms12/include $(ROOTCXXFLAGS)

ROOTLIBS     = $(shell $(ROOTSYS)/bin/root-config --libs) #-lRooFitCore

OBJS         = $(TMPDIR)/RooCMSShape.o \
	       $(TMPDIR)/RooErfXGaussian.o \
               $(TMPDIR)/TagAndProbe_dict.o

LIB          = libTagAndProbe.so


all: setup lib

setup:
	mkdir -p $(TMPDIR)
	mkdir -p $(LIBDIR)
	mkdir -p $(BINDIR)

lib: $(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(ROOTLIBS) -o $(LIBDIR)/$(LIB)

clean:
	rm -rf $(OBJS) $(LIBDIR)/$(LIB) \
	       $(TMPDIR)/TagAndProbe_dict.h $(TMPDIR)/TagAndProbe_dict.cc


################################################################################
# $(OBJS)
################################################################################

$(TMPDIR)/RooCMSShape.o: interface/RooCMSShape.h \
				    src/RooCMSShape.cc
	$(CXX) $(CXXFLAGS) -c src/RooCMSShape.cc \
	-o $(TMPDIR)/RooCMSShape.o 

$(TMPDIR)/RooErfXGaussian.o: interface/RooErfXGaussian.h \
				    src/RooErfXGaussian.cc
	$(CXX) $(CXXFLAGS) -c src/RooErfXGaussian.cc \
	-o $(TMPDIR)/RooErfXGaussian.o 

$(TMPDIR)/TagAndProbe_dict.o: $(TMPDIR)/TagAndProbe_dict.cc
	$(CXX) $(CXXFLAGS) -I$(TMPDIR) -I../.. -I/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms12/include -c $(TMPDIR)/TagAndProbe_dict.cc \
	-o $(TMPDIR)/TagAndProbe_dict.o

$(TMPDIR)/TagAndProbe_dict.cc: ../TagAndProbe/interface/RooCMSShape.h \
			       ../TagAndProbe/interface/RooErfXGaussian.h \
			   ../TagAndProbe/interface/PhysicsToolsTagAndProbeLinkDef.h
	rm -rf $(TMPDIR)/TagAndProbe_dict.h
	rm -rf $(TMPDIR)/TagAndProbe_dict.cc
	$(ROOTSYS)/bin/rootcint -f $(TMPDIR)/TagAndProbe_dict.cc \
	-c -I$(TMPDIR) -I../.. -I/afs/cern.ch/cms/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms12/include \
	../TagAndProbe/interface/RooCMSShape.h \
	../TagAndProbe/interface/RooErfXGaussian.h \
	../TagAndProbe/interface/PhysicsToolsTagAndProbeLinkDef.h