# Binaries
NAME = vrp
G_NAME = vrpgen

# Compiler
CC = g++

# Compiler flags
CXXFLAGS = -O3 -Wall

# COIN base directory
#COINDIR = /home/andrea/lnk/COIN/cbcDbg/
COINDIR = /home/andrea/lnk/COIN/coin-Cbc-2.7/
#COINDIR = <path>

# OSI libraries
OSIINC = -I$(COINDIR)include/coin
OSILIB = -L$(COINDIR)lib -lOsiCbc -lOsiClp -lCbc -lCbcSolver -lCgl \
                         -lClp -lOsi -lCoinUtils

# Libraries
LIB = $(OSILIB) -lm

# Set sources path
SRC_PATH = src
VPATH = $(SRC_PATH)

# Targets

default: $(NAME)
	
all: $(NAME) $(G_NAME)
	
	
# Clean

clean:
	rm -f *.o *~ $(NAME) $(G_NAME);
	cd src; rm -f *~; cd ..

# Generator build

GENOBJ = utils.o generator.o

$(G_NAME): $(GENOBJ)
	$(CC) -o $(G_NAME) $(GENOBJ) $(CXXFLAGS)

# Generator dependencies

generator.o: generator.cpp utils.hpp
	$(CC) -c $< -o $@ $(CXXFLAGS)


# Solver build

OBJ = VrpSolver.o VrpData.o VrpHeuristicSeparator.o utils.o main.o

$(NAME): $(OBJ)
	$(CC) -o $(NAME) $(OBJ) $(LIB) $(CXXFLAGS)

# Solver dependencies

main.o: main.cpp utils.hpp
	$(CC) -c $< -o $@ $(OSIINC) $(CXXFLAGS)

VrpSolver.o: VrpSolver.cpp utils.hpp VrpTypes.hpp
	$(CC) -c $< -o $@ $(OSIINC) $(CXXFLAGS)

utils.o: utils.cpp
	$(CC) -c $< -o $@ $(OSIINC) $(CXXFLAGS)

VrpData.o: VrpData.cpp VrpTypes.hpp
	$(CC) -c $< -o $@ $(OSIINC) $(CXXFLAGS)

VrpHeuristicSeparator.o: VrpHeuristicSeparator.cpp VrpTypes.hpp
	$(CC) -c $< -o $@ $(OSIINC) $(CXXFLAGS)

