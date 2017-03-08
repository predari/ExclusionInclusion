# move this makefile to the root of the project
DIR1:= /usr/lib
DIR2:= ./
INCDIR:= ./include
LIBDIR:= ./lib

program_NAME := myprogram
program_C_SRCS := $(wildcard *.c)
program_CXX_SRCS := $(wildcard *.cpp)
program_C_OBJS := ${program_C_SRCS:.c=.o}
program_CXX_OBJS := ${program_CXX_SRCS:.cpp=.o}
program_OBJS := $(program_C_OBJS) $(program_CXX_OBJS)
#program_INCLUDE_DIRS :=  # change to ./include once the project is ready
#program_LIBRARY_DIRS := 
#program_LIBRARY_DIRS := $(DIR1) $(DIR2)
program_LIBRARIES := lapacke m

CPPFLAGS += $(foreach includedir,$(program_INCLUDE_DIRS),-I$(includedir))
LDFLAGS += $(foreach librarydir,$(program_LIBRARY_DIRS),-L$(librarydir))
LDFLAGS += $(foreach library,$(program_LIBRARIES),-l$(library))

CXX := gcc
CXXFLAGS := -Wall
#LINK.cc := $(CXX) $(CPPFLAGS) $(LDFLAGS) 

.PHONY: all clean distclean

all: $(program_NAME)

$(program_NAME): $(program_OBJS)
	$(CXX) $(CXXFLAGS) $(program_OBJS) -o $(program_NAME) $(LDFLAGS)

clean:
	@- $(RM) $(program_NAME)
	@- $(RM) $(program_OBJS)

distclean: clean
