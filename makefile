# MAKEFILE for fRG_O2 Project

#--------------------------------------General------------------------------------------
CXX 	:= g++ # This is the main compiler
#CXX 	:= icc # This is the main compiler
SRCDIR 	:= src
HEADDIR := include
BUILDDIR := build
TARGETS := bin/gzero bin/anderson bin/hubbard2d

#--------------------------------------Sources and header files------------------------------------------
SRCEXT 	:= cpp
HEADEXT := h
SOURCES := $(shell ls -1 $(SRCDIR)/*.$(SRCEXT)) 
HEADERS := $(shell ls -1 $(HEADDIR)/*.$(HEADEXT)) 
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

#--------------------------------------Compiler settings------------------------------------------

CFLAGS += -std=c++11 -D OPENMP -fopenmp #			General compiler flags
CFLAGS += -D BOOST_DISABLE_ASSERTS #				Disable boost assert checks
DBFLAGS := -O2 -g #						Compiler flags for debugging
PROFFLAGS := -O2 -g #						Compiler flags for profiling
RUNFLAGS := -O3 #						Compiler flags for quick compile and run
OPTIMIZEFLAGS := -flto -march=native -O3 #			GCC Compiler flags for optimal speed
#OPTIMIZEFLAGS := -O3 -fp-model fast=2 -xHost # -no-prec-div #	Intel Compiler flags for optimal speed
H5LIB := $(shell h5c++ -show|cut -d ' ' -f2-) #			HDF5 libraries
LIB :=  -lmpi -lboost_serialization -lboost_mpi -lboost_timer $(H5LIB) -lpomerol -lmpi_cxx #		Library flags
INC := -I include -I/opt/pomerol/include -I/usr/include/hdf5/serial #			Additional include paths


#--------------------------------------Targets ------------------------------------------
run: 	CFLAGS += $(RUNFLAGS)
run: 	$(TARGETS)

debug: 	CFLAGS += $(DBFLAGS)
debug: 	$(TARGETS)

prof: 	CFLAGS += $(PROFFLAGS)
prof: 	$(TARGETS)

bin/anderson: prog/anderson.cpp
bin/hubbard2d: prog/hubbard2d.cpp
bin/gzero: prog/gzero.cpp
bin/dimer: prog/dimer.cpp

bin/%: $(OBJECTS)
	@mkdir -p bin
	@mkdir -p dat
	@echo " $(CXX) $^ -o $@ $(LIB) $(CFLAGS) $(INC)"; $(CXX) $^ -o $@ $(LIB) $(CFLAGS) $(INC)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) $(HEADERS)
	@mkdir -p $(BUILDDIR)
	@echo " $(CXX) $(CFLAGS) $(INC) -c -o $@ $<"; $(CXX) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGETS)

.PHONY: clean

