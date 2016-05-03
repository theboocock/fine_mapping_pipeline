
# -----------------------------------------------------------------
#   Makefile for GCTA 
#   
#   Supported platforms: Unix / Linux
# ---------------------------------------------------------------------

# Directory of the target
OUTPUT = gcta64

# Compiler
CXX = g++

# EIGEN library
#EIGEN_PATH = ../../../Lib/eigen

# Intel MKL library
#MKL_PATH = /opt/intel/mkl

# Compiler flags
CXXFLAGS = -w -O3 -m64 -static -fopenmp -I $(EIGEN_PATH) -DEIGEN_NO_DEBUG -I $(MKL_PATH)/include
LIB += -static -lz -Wl,--start-group  $(MKL_PATH)/lib/intel64/libmkl_intel_lp64.a $(MKL_PATH)/lib/intel64/libmkl_gnu_thread.a $(MKL_PATH)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

HDR += CommFunc.h \
	   cdflib.h \
	   dcdflib.h \
           gcta.h \
	   ipmpar.h \
           StatFunc.h \
           StrFunc.h \
           zfstream.h
SRC = bivar_reml.cpp \
           CommFunc.cpp \
           data.cpp \
	   dcdflib.cpp \
           est_hsq.cpp \
           grm.cpp \
           gwas_simu.cpp \
           ld.cpp \
           joint_meta.cpp \
           mlm_assoc.cpp \
	   mkl.cpp \
           option.cpp \
           popu_genet.cpp \
           raw_geno.cpp \
           StatFunc.cpp \
           StrFunc.cpp \
           zfstream.cpp
	   
OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT) 

$(OUTPUT) :
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean: 
	rm -f *.o
