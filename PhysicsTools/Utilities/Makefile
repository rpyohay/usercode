################################################################################
#
# PhysicsTools/Utilities Makefile (for standalone use outside CMSSW/SCRAM)
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
CXXFLAGS     = -O3 -Wall -fPIC -DSTANDALONE -I../.. $(ROOTCXXFLAGS)

ROOTLIBS     = $(shell $(ROOTSYS)/bin/root-config --libs)

OBJS         = $(TMPDIR)/LumiReweightingStandAlone.o \
               $(TMPDIR)/Utilities_dict.o

LIB          = libUtilities.so


all: setup lib

setup:
	mkdir -p $(TMPDIR)
	mkdir -p $(LIBDIR)
	mkdir -p $(BINDIR)

lib: $(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(ROOTLIBS) -o $(LIBDIR)/$(LIB)

clean:
	rm -rf $(OBJS) $(LIBDIR)/$(LIB) \
	       $(TMPDIR)/Utilities_dict.h $(TMPDIR)/Utilities_dict.cc


################################################################################
# $(OBJS)
################################################################################

$(TMPDIR)/LumiReweightingStandAlone.o: interface/LumiReweightingStandAlone.h \
				    src/LumiReweightingStandAlone.cc
	$(CXX) $(CXXFLAGS) -c src/LumiReweightingStandAlone.cc \
	-o $(TMPDIR)/LumiReweightingStandAlone.o 

$(TMPDIR)/Utilities_dict.o: $(TMPDIR)/Utilities_dict.cc
	$(CXX) $(CXXFLAGS) -I$(TMPDIR) -I../.. -c $(TMPDIR)/Utilities_dict.cc \
	-o $(TMPDIR)/Utilities_dict.o

$(TMPDIR)/Utilities_dict.cc: ../Utilities/interface/LumiReweightingStandAlone.h \
			   ../Utilities/interface/Linkdef.h
	rm -rf $(TMPDIR)/Utilities_dict.h
	rm -rf $(TMPDIR)/Utilities_dict.cc
	$(ROOTSYS)/bin/rootcint -f $(TMPDIR)/Utilities_dict.cc \
	-c -I$(TMPDIR) -I../.. \
	../Utilities/interface/LumiReweightingStandAlone.h \
	../Utilities/interface/Linkdef.h