CC=gcc
CXX=g++
RM=rm -f
CPPFLAGS=-Wall -std=gnu++17 -O2 -isystem /usr/include/eigen3
LDLIBS=-lglpk

SRCS=LPSolver.cpp Problem.cpp find_constraint.cpp
OBJS=$(patsubst %.cpp,bin/%.o,$(SRCS))

all: bin/clqo

bin/clqo: bin/test.o lib
	$(CXX) -o bin/clqo $(OBJS) bin/test.o $(LDLIBS) 

lib: $(OBJS)

bin/%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c $< -o $@ 

clean:
	$(RM) bin/*.o

distclean: clean
	$(RM) tool