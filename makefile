CC=gcc
CXX=g++
RM=rm -f
CPPFLAGS=-Wall -std=c++0x -O2 -I /usr/include/eigen3
LDLIBS=-lglpk

SRCS=LPSolver.cpp Problem.cpp 
OBJS=$(patsubst %.cpp,bin/%.o,$(SRCS))

all: bin/clqo

bin/clqo: bin/test.o lib
	$(CXX) -o bin/clqo $(OBJS) test.o $(LDLIBS) 

lib: bin/LPSolver.o bin/Problem.o

bin/%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c $< -o $@ 

clean:
	$(RM) bin/*.o

distclean: clean
	$(RM) tool