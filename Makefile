
CXXFLAGS = -O3 -msse3 -g $(INCLUDE)
LDFLAGS = -lboost_program_options

all: simulation

simulation: simulation.cpp

clean:
	rm -Rf simulation *.o
