# Makefile for excusionInclusion tool.
# Depend from lapack which needs blas
# add in LDFLAGS the path to lapack
# include the path where lapacke is
# link with lapacke (-llapacke)


# lapack directory
LAPACK_DIR:= /usr/local/opt/lapack/lib
# not used 
# INCDIR:= ./include
# LIBDIR:= ./lib

program_NAME := myprogram
program_C_SRCS := $(wildcard *.c)
program_CXX_SRCS := $(wildcard *.cpp)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
program_INCLUDE_DIRS := /usr/local/Cellar/lapack/3.8.0/include
program_LIBRARY_DIRS := 
program_LIBRARY_DIRS := $(LAPACK_DIR)
program_LIBRARIES := lapacke m

CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))

CXX := gcc
CXXFLAGS := -Wall -std=c99
#LINK.cc := $(CXX) $(CPPFLAGS) $(LDFLAGS) 

.PHONY: all clean distclean

all: $(program_NAME)

$(program_NAME): $(program_OBJS)
	$(CXX) $(CXXFLAGS) $(program_OBJS) -o $(program_NAME) $(LDFLAGS)

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)

distclean: clean
