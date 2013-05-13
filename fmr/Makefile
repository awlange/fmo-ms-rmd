# ----------------------------------------------- #
# ** Makefile for FMR code **
# ----------------------------------------------- #


# ----------------------------------------------- #
# ** User defined variables **
# ----------------------------------------------- #
# MPI C++ compiler
CPP = mpicxx

# Debugging flags
#DEBUG = -g -DFMR_DEBUG
DEBUG = -g 

# Compiler flags
CFLAGS = -O3 $(DEBUG)

# Libraries, must link against math and LAPACK libraries 
LIBDIR = -L/usr/lib -L/usr/lib/lapack
LIB = -lm -llapack

# ----------------------------------------------- #
# Compilation instructions
# ----------------------------------------------- #

SRC = $(wildcard *.cpp)
INC = $(wildcard *.h)
OBJ = $(SRC:.cpp=.o)

install: $(OBJ)
	$(CPP) $(OBJ) $(LIB) -o fmr.exe

clean: 
	rm -f *.o fmr.exe

%.o:%.cpp
	$(CPP) $(CFLAGS) -c $<
