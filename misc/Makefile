# Makefile for excusionInclusion tool.
# Depends from lapack which needs blas
# add in LDFLAGS the path to lapack
# include the path where lapacke is
# link with lapacke (-llapacke)


# lapack directory
LAPACK_DIR:= /usr/local/opt/lapack/lib
#LAPACK_DIR:=/usr/local/Cellar/lapack/3.8.0/lib
# not used 
# INCDIR:= ./include
# LIBDIR:= ./lib

program1_NAME := ts
program2_NAME := psc
#program_C_SRCS := $(wildcard *.c)
program_C_SRCS := aux.c
program_CXX_SRCS := $(wildcard *.cpp)
program_C_OBJS := ${program_C_SRCS:.c=.o}	
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
program1 := test.c
program2 := pseudoCompare.c
#program_INCLUDE_DIRS := /usr/local/Cellar/lapack/3.8.0/include
program_INCLUDE_DIRS := /usr/local/opt/lapack/include
program_LIBRARY_DIRS := 
program_LIBRARY_DIRS := $(LAPACK_DIR)
program_LIBRARIES := lapacke m #blas lapack gfortran

CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))

# for my mac 
#CXX := gcc-8
# for my linux
CXX := gcc
CXXFLAGS := -Wall -std=c99
LINK.cc := $(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) 

.PHONY: all clean distclean

all: $(program1_NAME) $(program2_NAME)

%.o: %.c aux.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $<

$(program1_NAME): $(program_OBJS)
	$(CXX) $(CXXFLAGS) $(program1) $(program_C_OBJS) -o $(program1_NAME) $(LDFLAGS)

$(program2_NAME): $(program_OBJS)
	$(CXX) $(CXXFLAGS) $(program2) $(program_C_OBJS) -o $(program2_NAME) $(LDFLAGS)

# all: program1 program2
# program1: program1.c
#     gcc -o program1 program1.c
# program2: program2.c
#     gcc -o program2 program2.c


clean:
	@- $(RM) $(program1_NAME)
	@- $(RM) $(program2_NAME)
	@- $(RM) $(program_OBJS)

distclean: clean




