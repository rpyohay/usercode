################################################################################
#
# GMSBTools/Filters Makefile (for standalone use outside CMSSW/SCRAM)
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

OBJS         = $(TMPDIR)/Categorizer.o \
               $(TMPDIR)/Filters_dict.o

LIB          = libFilters.so


all: setup lib

setup:
	mkdir -p $(TMPDIR)
	mkdir -p $(LIBDIR)
	mkdir -p $(BINDIR)

lib: $(OBJS)
	$(CXX) $(CXXFLAGS) -shared $(OBJS) $(ROOTLIBS) -o $(LIBDIR)/$(LIB)

clean:
	rm -rf $(OBJS) $(LIBDIR)/$(LIB) CondFormats \
	rm -rf $(OBJS) $(LIBDIR)/$(LIB) \
	       $(TMPDIR)/Filters_dict.h $(TMPDIR)/Filters_dict.cc


################################################################################
# $(OBJS)
################################################################################

$(TMPDIR)/Categorizer.o: interface/Categorizer.h \
			 interface/Typedefs.h \
				    src/Categorizer.cc
	$(CXX) $(CXXFLAGS) -c src/Categorizer.cc \
	-o $(TMPDIR)/Categorizer.o 

$(TMPDIR)/Filters_dict.o: $(TMPDIR)/Filters_dict.cc
	$(CXX) $(CXXFLAGS) -I$(TMPDIR) -I../.. -I$(INC) -c $(TMPDIR)/Filters_dict.cc \
	-o $(TMPDIR)/Filters_dict.o

$(TMPDIR)/Filters_dict.cc: ../Filters/interface/Categorizer.h \
			   ../Filters/interface/Linkdef.h
	rm -rf $(TMPDIR)/Filters_dict.h
	rm -rf $(TMPDIR)/Filters_dict.cc
	$(ROOTSYS)/bin/rootcint -f $(TMPDIR)/Filters_dict.cc \
	-c -I$(TMPDIR) -I../.. -I$(INC) \
	../Filters/interface/Categorizer.h \
	../Filters/interface/Typedefs.h \
	../Filters/interface/Linkdef.h
