# MAKEFILE for fRG_O2 Project

#--------------------------------------General------------------------------------------
CXX 	:= g++ # This is the main compiler
#CXX 	:= icc # This is the main compiler
SRCDIR 	:= src
HEADDIR := include
BUILDDIR := build
TARGET 	:= bin/run

#--------------------------------------Sources and header files------------------------------------------
SRCEXT 	:= cpp
HEADEXT := h
SOURCES := $(shell ls -1 $(SRCDIR)/*.$(SRCEXT)) 
HEADERS := $(shell ls -1 $(HEADDIR)/*.$(HEADEXT)) 
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

#--------------------------------------Compiler settings------------------------------------------

CFLAGS += -std=c++11 # 								General compiler flags
OPTIMIZEFLAGS := -march=native -O3 #-flto # 					Compiler flags for optimal speed
DBFLAGS := -O2 -g # 								Compiler flags for debugging
#OPTIMIZEFLAGS := -xHOST -fp-model fast -lpthread -O3 -openmp -DOPENMP # -opt-prefetch 0 # Compiler flags for optimal speed
LIB :=  -lmpi -lboost_serialization -lboost_mpi -lpomerol -lmpi_cxx #		Library flags
INC := -I include #-I/opt/pomerol/include #						Additional include paths

#--------------------------------------Targets ------------------------------------------

all: 	CFLAGS += $(OPTIMIZEFLAGS)
all: 	bin/anderson bin/gzero

debug: 	CFLAGS += $(DBFLAGS)
debug: 	bin/anderson bin/gzero

bin/%: $(SRCDIR)/%.$(SRCEXT) $(HEADERS)
	@mkdir -p $(BUILDDIR)
	@mkdir -p bin
	@mkdir -p log
	@mkdir -p dat
	@echo " $(CXX) -o $@ $< $(CFLAGS) $(INC) $(LIB)"; $(CXX) -o $@ $< $(CFLAGS) $(INC) $(LIB)

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean

