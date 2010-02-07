
PLATFORM = $(shell uname)

ifeq ($(PLATFORM), Darwin)
	INCLUDE = -I/software/local/include
else
	INCLUDE = 
endif

CXXFLAGS = -O3 -msse3 -g $(INCLUDE)

all: simulation

simulation: simulation.cpp

clean:
	rm -Rf simulation *.o *.dSYM
